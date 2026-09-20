import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from run import compare

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

if __name__=='__main__':unittest.main()
