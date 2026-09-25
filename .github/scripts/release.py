"""Resolve immutable release inputs and package the tested plugin (stdlib only)."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
from urllib.parse import quote
import xml.etree.ElementTree as ET


def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args], text=True).strip()


def version(source_text):
    match = re.search(r'project\(neo_fft\s+VERSION\s+(\d+\.\d+\.\d+)', source_text)
    if not match:
        raise ValueError('Cannot read neo_fft project version')
    return match[1]


def resolve():
    kind = os.environ['GITHUB_REF_TYPE']
    ref = os.environ['GITHUB_REF_NAME']
    commit = git('.', 'rev-parse', '--verify', f"{os.environ['GITHUB_SHA']}^{{commit}}")
    if git('.', 'rev-parse', 'HEAD') != commit:
        raise ValueError('Source checkout differs from the workflow run commit')
    if kind == 'tag':
        if not re.fullmatch(r'v?\d+\.\d+\.\d+(?:-[0-9A-Za-z]+(?:[.-][0-9A-Za-z]+)*)?', ref):
            raise ValueError('Publishing requires a version tag')
        actual = version(git('.', 'show', f'{commit}:CMakeLists.txt'))
        if ref.removeprefix('v').split('-', 1)[0] != actual:
            raise ValueError(f'Tag {ref} does not match source version {actual}')
        values = dict(commit=commit, kind=kind, tag=ref, label=ref)
    elif kind == 'branch':
        values = dict(commit=commit, kind=kind, tag='', label=f'git-{commit[:12]}')
    else:
        raise ValueError('Select a branch or tag when dispatching the workflow')
    with open(os.environ['GITHUB_OUTPUT'], 'a', encoding='utf-8') as output:
        for key, value in values.items():
            output.write(f'{key}={value}\n')
    with open(os.environ['GITHUB_STEP_SUMMARY'], 'a', encoding='utf-8') as output:
        output.write(f"Source: `{commit}`\n\nMode: **{values['kind']}**; "
                     f"{'publish after all six builds pass' if kind == 'tag' else 'artifacts only; no GitHub Release'}.\n")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def package_docs(source, stage, commit):
    roots = [source / 'docs' / section for section in ('api', 'knowledge')]
    base = f"{os.environ['GITHUB_SERVER_URL']}/{os.environ['GITHUB_REPOSITORY']}/blob/{commit}/"
    for root in roots:
        dest = stage / root.relative_to(source)
        shutil.copytree(root, dest)
        for page in dest.rglob('*.md'):
            original = source / page.relative_to(stage)

            def link(match):
                target = match.group(1)
                if '://' in target or target.startswith('#'):
                    return match.group(0)
                path, sep, anchor = target.partition('#')
                resolved = (original.parent / path).resolve()
                if any(resolved.is_relative_to(r) for r in roots):
                    return match.group(0)
                if not resolved.is_file():
                    raise ValueError(f'Missing documentation target: {target}')
                # Source references resolve at this exact commit, not a moving branch.
                url = base + quote(resolved.relative_to(source).as_posix())
                return '](' + url + (sep + anchor if sep else '') + ')'

            page.write_text(re.sub(r'\]\(([^)]+)\)', link, page.read_text(encoding='utf-8')),
                            encoding='utf-8')


def package(args):
    source, build, output = (Path(p).resolve() for p in (args.source, args.build, args.output))
    commit = os.environ['RELEASE_COMMIT']
    if git(source, 'rev-parse', 'HEAD') != commit:
        raise ValueError('Source checkout differs from resolved commit')
    tests_path = build / 'tests.xml'
    tests = list(ET.parse(tests_path).getroot().iter('testcase'))
    if not tests or any(t.find('failure') is not None or t.find('error') is not None
                        or t.find('skipped') is not None or t.get('status') == 'notrun' for t in tests):
        raise ValueError('All configured tests must pass without skips before packaging')
    host_tests = os.environ['RELEASE_HOST_TESTS'] == 'ON'
    if host_tests and not any(t.get('name', '').startswith('vs.') for t in tests):
        raise ValueError('VapourSynth host tests were requested but are absent from test evidence')
    candidates = [build / 'Release/neo-fft.dll', build / 'neo-fft.so', build / 'neo-fft.dylib']
    binaries = [p for p in candidates if p.is_file()]
    if len(binaries) != 1:
        raise ValueError(f'Expected one release plugin, found {binaries}')
    name = f"neo-fft-{os.environ['RELEASE_LABEL']}-{os.environ['RELEASE_PLATFORM']}"
    stage = build / 'package' / name
    stage.mkdir(parents=True, exist_ok=False)
    binary = binaries[0]
    shutil.copy2(binary, stage / binary.name)
    shutil.copy2(tests_path, stage / 'tests.xml')
    package_docs(source, stage, commit)
    recipes = stage / 'build-recipes'
    (recipes / 'cmake').mkdir(parents=True)
    shutil.copy2(source / 'CMakeLists.txt', recipes / 'CMakeLists.txt')
    shutil.copy2(source / 'cmake/Dependencies.cmake', recipes / 'cmake/Dependencies.cmake')
    dependencies = {}
    for dep in ('dualsynth2', 'highway', 'pocketfft'):
        dep_root = build / '_deps' / f'{dep}-src'
        licenses = sorted(p for p in dep_root.iterdir()
                          if p.is_file() and p.name.upper().startswith(('LICENSE', 'COPYING')))
        if any(not p.is_file() for p in licenses):
            raise ValueError(f'Missing dependency notice file: {dep}')
        if not licenses:
            raise ValueError(f'Missing dependency notices: {dep}')
        dest = stage / 'licenses' / dep
        dest.mkdir(parents=True)
        for path in licenses:
            shutil.copy2(path, dest / path.name)
        dependencies[dep] = {'notices': {p.name: sha(p) for p in licenses}}
        if (dep_root / '.git').exists():
            dependencies[dep]['commit'] = git(dep_root, 'rev-parse', 'HEAD')
    for path in source.iterdir():
        if path.is_file() and path.name.upper().startswith(('LICENSE', 'COPYING')):
            shutil.copy2(path, stage / path.name)
    manifest = {
        'commit': commit, 'version': version((source / 'CMakeLists.txt').read_text()),
        'tag': os.environ['RELEASE_TAG'] or None, 'mode': os.environ['RELEASE_KIND'],
        'platform': os.environ['RELEASE_PLATFORM'], 'runner': platform.platform(),
        'workflow_commit': os.environ['RELEASE_WORKFLOW_COMMIT'],
        'run_url': f"{os.environ['GITHUB_SERVER_URL']}/{os.environ['GITHUB_REPOSITORY']}/actions/runs/{os.environ['GITHUB_RUN_ID']}",
        'build_type': 'Release', 'simd': True, 'tests_passed': len(tests),
        'vapoursynth_host_tests': host_tests, 'avisynth_host_tests': False,
        'dependencies': dependencies,
        'plugin': {'file': binary.name, 'sha256': sha(binary)},
    }
    (stage / 'build-info.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    (stage / 'README.txt').write_text(
        f"neo-fft {manifest['version']} ({commit})\nPlatform: {manifest['platform']}\n\n"
        "Install the plugin in your matching VapourSynth plugin directory, or load it explicitly.\n"
        "For AviSynth interface 11 or later, load this same plugin with LoadPlugin; functions use neo_fft_.\n"
        "This archive does not include either host runtime. AviSynth host tests are not run by this workflow.\n"
        "Linux binaries are built on Ubuntu 24.04; other distributions require compatible runtime libraries.\n"
        "macOS binaries are built on macOS 15 and are not signed or notarized.\n\n"
        f"Spec/SIMD tests passed: {len(tests)}. VapourSynth host tests: "
        f"{'passed' if host_tests else 'not configured on this platform'}.\n"
        "See docs/api/ and docs/knowledge/ for the bilingual API reference and knowledge base.\n"
        "See build-info.json, tests.xml, build-recipes/, and licenses/ for build evidence and dependency notices.\n",
        encoding='utf-8')
    output.mkdir(parents=True, exist_ok=True)
    archive = Path(shutil.make_archive(str(output / name), 'zip', stage.parent, stage.name))
    (output / f'{name}.sha256').write_text(f'{sha(archive)}  {archive.name}\n', encoding='ascii')
    print(archive)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('resolve')
    packaging = commands.add_parser('package')
    packaging.add_argument('--source', required=True)
    packaging.add_argument('--build', required=True)
    packaging.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.command == 'resolve':
        resolve()
    else:
        package(args)
