import contextlib
import io
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from run import compare, reuse_candidate

class Protocol(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.cases=[]
        self.budgets={'limits':{}}
        for algorithm in ('FFT3D','DFTTest'):
            case=dict(id=algorithm,algorithm=algorithm,bits=8,params={})
            self.cases.append(case)
            for label in ('old1','new1','old0','new0'):
                path=self.root/algorithm/label;path.mkdir(parents=True)
                record=dict(case=case,input_sha256='fixture-hash',status='ok')
                manifest=dict(kernel_info={'target':'unit-test-vector-target'},cases=[record])
                (path/'manifest.json').write_text(json.dumps(manifest))
                np.savez(path/(algorithm+'.npz'),**{'0-0':np.zeros((2,2),dtype=np.uint8)})
            for pair in 'ABCD': self.budgets['limits'][f'{algorithm}/8/pocketfft/{pair}']={'atol':0}
    def change(self,label,edit):
        path=self.root/'FFT3D'/label/'manifest.json'
        value=json.loads(path.read_text());edit(value);path.write_text(json.dumps(value))
    def compare(self):
        with contextlib.redirect_stdout(io.StringIO()): return compare(self.root,self.cases,self.budgets)
    def test_equal(self): self.assertFalse(self.compare())
    def test_stale_output_cannot_mask_capture_error(self):
        self.change('new1',lambda v:v['cases'][0].update(status='error'))
        self.assertTrue(self.compare())
    def test_input_hash_mismatch(self):
        self.change('new1',lambda v:v['cases'][0].update(input_sha256='different'))
        with self.assertRaisesRegex(RuntimeError,'inputs'): self.compare()
    def test_scalar_is_not_simd_evidence(self):
        self.change('new0',lambda v:v.update(kernel_info={'target':'SCALAR'}))
        with self.assertRaisesRegex(RuntimeError,'scalar target'): self.compare()
    def test_budget_applies_to_each_sample(self):
        np.savez(self.root/'FFT3D/new1/FFT3D.npz',**{'0-0':np.ones((2,2),dtype=np.uint8)})
        self.assertTrue(self.compare())
    def test_public_reference_requires_autoload_and_fftw_evidence(self):
        with self.assertRaisesRegex(RuntimeError,'autoload/FFTW'):
            compare(self.root,self.cases,None,'fftw-pocketfft')
    def test_public_reference_uses_separate_backend_budget(self):
        for algorithm in ('FFT3D','DFTTest'):
            for label in ('old1','old0'):
                path=self.root/algorithm/label/'manifest.json'
                value=json.loads(path.read_text())
                value.update(autoload=True,fft_backend='fftw',loaded_fft_libraries=[{'path':'fixture-fftw'}])
                path.write_text(json.dumps(value))
        with self.assertRaisesRegex(RuntimeError,'missing frozen budget.*fftw-pocketfft'):
            compare(self.root,self.cases,self.budgets,'fftw-pocketfft')
        limits={k.replace('/pocketfft/','/fftw-pocketfft/'):v for k,v in self.budgets['limits'].items()}
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(compare(self.root,self.cases,{'limits':limits},'fftw-pocketfft'))
    def reuse_inputs(self):
        source=self.root/'FFT3D/new1'; plugin=self.root/'candidate.dll'; plugin.write_bytes(b'candidate fixture')
        path=source/'manifest.json'; value=json.loads(path.read_text())
        value.update(sha256=hashlib.sha256(plugin.read_bytes()).hexdigest(),seed=17,opt=1)
        value['cases'][0]['output_hashes']=[hashlib.sha256(np.zeros((2,2),dtype=np.uint8).tobytes()).hexdigest()]
        path.write_text(json.dumps(value))
        return source,self.root/'reuse',plugin,17,1,[self.cases[0]]
    def test_reuse_rejects_different_binary(self):
        args=self.reuse_inputs(); args[2].write_bytes(b'changed candidate')
        with self.assertRaisesRegex(RuntimeError,'binary'): reuse_candidate(*args)
    def test_reuse_rejects_changed_output(self):
        args=self.reuse_inputs()
        np.savez(args[0]/'FFT3D.npz',**{'0-0':np.ones((2,2),dtype=np.uint8)})
        with self.assertRaisesRegex(RuntimeError,'output hash'): reuse_candidate(*args)
    def test_reuse_records_provenance(self):
        args=self.reuse_inputs(); reuse_candidate(*args)
        value=json.loads((args[1]/'manifest.json').read_text())
        self.assertEqual(value['reused_from'],str(args[0].resolve()))
        self.assertEqual(value['original_manifest_sha256'],hashlib.sha256((args[0]/'manifest.json').read_bytes()).hexdigest())

if __name__=='__main__':unittest.main()
