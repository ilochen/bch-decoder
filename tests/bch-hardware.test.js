import test from 'node:test';import assert from 'node:assert/strict';
import {codes,field,syndromes,decode} from '../assets/models.js';
import {experiment,traceWord} from '../assets/bch-hardware.js';
for(const n of [63,255,1023]){
 test(`BCH ${n}: byte-wide Horner updates match independent polynomial evaluation`,()=>{
  const e=experiment(n),tr=e.candidates[0];assert.deepEqual(tr.S,syndromes(e.received));
  assert.deepEqual([...tr.roots].sort((a,b)=>a-b),[3,17]);assert.ok(e.recovered);
 });
 test(`BCH ${n}: scaled RTL-style BM and eight-way Chien preserve roots`,()=>{
  const f=field(n);
  for(let count=1;count<=codes[n].t;count++){
   const errors=Array.from({length:count},(_,i)=>7+i*13),e=experiment(n,'two','hard',errors),tr=e.candidates[0],reference=decode(e.received);
   assert.ok(tr.valid);assert.deepEqual([...tr.roots].sort((a,b)=>a-b),reference.roots);
   for(const frame of tr.frames.filter(x=>x.phase==='chien'))for(const lane of frame.lanes)if(lane.position>=0){
    const direct=tr.bm.sigma.reduce((sum,c,j)=>sum^f.mul(c,f.pow(-lane.position*j)),0);assert.equal(lane.value,direct);
   }
   assert.equal(tr.bm.trace.length,4*codes[n].t);
   for(const frame of tr.bm.trace.filter(x=>x.phase==='update'))assert.equal(frame.nextSigma[0],f.mul(frame.gamma,frame.sigma[0]));
  }
 });
 test(`BCH ${n}: clean word bypasses BM/search and soft decoding rescues t+1 errors`,()=>{
  const clean=experiment(n,'clean').candidates[0];assert.ok(clean.clean);assert.equal(clean.frames.length,Math.ceil(n/8)+1);
  const soft=experiment(n,'soft','soft');assert.ok(soft.recovered);assert.equal(soft.best.mask,1);
  const S=soft.candidates.map(c=>c.S);for(let i=0;i<S[0].length;i++)assert.equal(S[3][i],S[0][i]^S[1][i]^S[2][i]);
 });
}
test('zero position and last partial Chien group are not lost or duplicated',()=>{
 const bits=Array(63).fill(0);bits[0]=1;bits[62]=1;const tr=traceWord(bits);assert.deepEqual(tr.roots,[62,0]);assert.ok(tr.valid);
});
