"""Exercise manual branch/tag release selection in an isolated Git repository."""

import os
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile


SCRIPT = Path(__file__).with_name('release.py').resolve()


class ResolveTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.git('init', '-b', 'master')
        self.git('config', 'user.name', 'Release test')
        self.git('config', 'user.email', 'test@example.invalid')
        self.git('config', 'commit.gpgsign', 'false')
        self.git('config', 'tag.gpgsign', 'false')
        self.git('config', 'core.autocrlf', 'false')
        (self.root / 'CMakeLists.txt').write_text('project(neo_fft VERSION 0.9.0)\n')
        self.git('add', 'CMakeLists.txt')
        self.git('commit', '-m', 'base')
        self.base = self.git('rev-parse', 'HEAD')
        self.git('tag', '0.9.0')

    def git(self, *args):
        return subprocess.check_output(['git', '-C', str(self.root), *args],
                                       text=True, stderr=subprocess.PIPE).strip()

    def resolve(self, kind, ref, sha=None):
        output = self.root / 'output'
        output.write_text('')
        env = dict(os.environ, GITHUB_REF_TYPE=kind, GITHUB_REF_NAME=ref,
                   GITHUB_SHA=sha or self.base,
                   GITHUB_OUTPUT=str(output), GITHUB_STEP_SUMMARY=str(self.root / 'summary'))
        result = subprocess.run([sys.executable, str(SCRIPT), 'resolve'], cwd=self.root,
                                env=env, capture_output=True, text=True)
        return result, dict(line.split('=', 1) for line in output.read_text().splitlines())

    def test_branches_only_produce_artifacts(self):
        for ref in ('master', 'feature/test', '0.9.0'):
            with self.subTest(ref=ref):
                result, output = self.resolve('branch', ref)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(output, dict(commit=self.base, kind='branch', tag='',
                                              label=f'git-{self.base[:12]}'))

    def test_version_tags_publish(self):
        for ref in ('0.9.0', 'v0.9.0', '0.9.0-rc.1'):
            with self.subTest(ref=ref):
                result, output = self.resolve('tag', ref)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(output, dict(commit=self.base, kind='tag', tag=ref, label=ref))

    def test_annotated_tag_resolves_to_commit(self):
        self.git('tag', '-a', 'v0.9.0', '-m', 'release')
        result, output = self.resolve('tag', 'v0.9.0', self.git('rev-parse', 'v0.9.0'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output['commit'], self.base)

    def test_invalid_or_mismatched_tag_fails(self):
        for ref in ('nightly', '1.0.0', '0.9.0\nother=value'):
            with self.subTest(ref=ref):
                result, output = self.resolve('tag', ref)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(output, {})

    def test_checkout_must_match_run_commit(self):
        self.git('commit', '--allow-empty', '-m', 'later commit')
        result, output = self.resolve('branch', 'master')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('differs from the workflow run commit', result.stderr)
        self.assertEqual(output, {})

    def test_branch_movement_does_not_change_run_commit(self):
        self.git('commit', '--allow-empty', '-m', 'branch advances')
        self.git('checkout', '--detach', self.base)
        result, output = self.resolve('branch', 'master')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output['commit'], self.base)


class PackageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'source'
        self.build = self.root / 'build'
        self.output = self.root / 'dist'
        self.source.mkdir()
        self.build.mkdir()
        for section in ('api', 'knowledge'):
            for language in ('en', 'zh-CN'):
                folder = self.source / 'docs' / section / language
                folder.mkdir(parents=True)
                (folder / 'README.md').write_text('[Source](../../../src/example.cpp)\n', encoding='utf-8')
        (self.source / 'src').mkdir()
        (self.source / 'src/example.cpp').write_text('// source reference\n')
        (self.source / 'cmake').mkdir()
        (self.source / 'cmake/Dependencies.cmake').write_text('# pinned dependencies\n')
        (self.source / 'CMakeLists.txt').write_text('project(neo_fft VERSION 0.9.0)\n')
        self.git('init')
        self.git('config', 'user.name', 'Release test')
        self.git('config', 'user.email', 'test@example.invalid')
        self.git('config', 'commit.gpgsign', 'false')
        self.git('config', 'core.autocrlf', 'false')
        self.git('add', '.')
        self.git('commit', '-m', 'fixture')
        self.commit = self.git('rev-parse', 'HEAD')
        for dep in ('dualsynth2', 'highway', 'pocketfft'):
            folder = self.build / '_deps' / f'{dep}-src'
            folder.mkdir(parents=True)
            (folder / 'LICENSE').write_text(f'{dep} notice\n')
        (self.build / 'tests.xml').write_text('<testsuite><testcase name="core.fixture" status="run"/></testsuite>')
        # A synthetic payload tests packaging only, not executable compatibility.
        self.binary = self.build / 'neo-fft.so'
        self.binary.write_bytes(b'packaging fixture')
        self.env = dict(os.environ, RELEASE_COMMIT=self.commit, RELEASE_LABEL=f'git-{self.commit[:12]}',
                        RELEASE_PLATFORM='linux-x64', RELEASE_HOST_TESTS='OFF', RELEASE_TAG='',
                        RELEASE_KIND='branch', RELEASE_WORKFLOW_COMMIT=self.commit,
                        GITHUB_SERVER_URL='https://github.com', GITHUB_REPOSITORY='example/neo-fft',
                        GITHUB_RUN_ID='123')

    def git(self, *args):
        return subprocess.check_output(['git', '-C', str(self.source), *args],
                                       text=True, stderr=subprocess.PIPE).strip()

    def package(self):
        return subprocess.run([sys.executable, str(SCRIPT), 'package', '--source', str(self.source),
                               '--build', str(self.build), '--output', str(self.output)],
                              env=self.env, capture_output=True, text=True)

    def test_package_contents_and_checksums(self):
        result = self.package()
        self.assertEqual(result.returncode, 0, result.stderr)
        archive, = self.output.glob('*.zip')
        checksum, = self.output.glob('*.sha256')
        digest, name = checksum.read_text().split()
        self.assertEqual(name, archive.name)
        self.assertEqual(digest, hashlib.sha256(archive.read_bytes()).hexdigest())
        prefix = archive.stem + '/'
        with zipfile.ZipFile(archive) as z:
            manifest = json.loads(z.read(prefix + 'build-info.json'))
            self.assertEqual(manifest['commit'], self.commit)
            self.assertEqual(manifest['version'], '0.9.0')
            self.assertFalse(manifest['vapoursynth_host_tests'])
            self.assertFalse(manifest['avisynth_host_tests'])
            self.assertEqual(manifest['tests_passed'], 1)
            self.assertEqual(manifest['plugin']['sha256'], hashlib.sha256(self.binary.read_bytes()).hexdigest())
            self.assertEqual(set(manifest['dependencies']), {'dualsynth2', 'highway', 'pocketfft'})
            for dep, record in manifest['dependencies'].items():
                data = z.read(prefix + f'licenses/{dep}/LICENSE')
                self.assertEqual(record['notices']['LICENSE'], hashlib.sha256(data).hexdigest())
            for section in ('api', 'knowledge'):
                for language in ('en', 'zh-CN'):
                    doc = z.read(prefix + f'docs/{section}/{language}/README.md').decode('utf-8')
                    self.assertIn(f'https://github.com/example/neo-fft/blob/{self.commit}/src/example.cpp', doc)
            self.assertIn(prefix + 'tests.xml', z.namelist())
            self.assertIn(prefix + 'build-recipes/CMakeLists.txt', z.namelist())
            self.assertIn(prefix + 'build-recipes/cmake/Dependencies.cmake', z.namelist())

    def test_reject_failed_skipped_or_empty_evidence(self):
        for body in ('', '<testcase><failure/></testcase>', '<testcase><error/></testcase>',
                     '<testcase><skipped/></testcase>', '<testcase status="notrun"/>'):
            with self.subTest(body=body):
                (self.build / 'tests.xml').write_text(f'<testsuite>{body}</testsuite>')
                result = self.package()
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('All configured tests must pass', result.stderr)

    def test_reject_missing_or_ambiguous_binary(self):
        self.binary.unlink()
        self.assertNotEqual(self.package().returncode, 0)
        self.binary.write_bytes(b'fixture')
        (self.build / 'neo-fft.dylib').write_bytes(b'other')
        self.assertNotEqual(self.package().returncode, 0)

    def test_reject_wrong_source_commit(self):
        self.env['RELEASE_COMMIT'] = '0' * 40
        result = self.package()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Source checkout differs', result.stderr)

    def test_reject_missing_dependency_notice(self):
        (self.build / '_deps/pocketfft-src/LICENSE').unlink()
        result = self.package()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Missing dependency notices', result.stderr)

    def test_reject_host_coverage_without_host_tests(self):
        self.env['RELEASE_HOST_TESTS'] = 'ON'
        result = self.package()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('absent from test evidence', result.stderr)


if __name__ == '__main__':
    unittest.main()
