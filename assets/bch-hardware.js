import {codes,field,encode,syndromes} from './models.js';
// Functional, register-level teaching trace. Control fill/drain and output handshakes are omitted.
export function inversionless(S,n){
 const f=field(n),trace=[];let sigma=[1,0,0,0,0],B=[1,0,0,0,0],gamma=1,L=0;
 for(let r=0;r<S.length;r++){
  const window=Array.from({length:5},(_,j)=>S[r-j]||0),products=sigma.map((v,j)=>f.mul(v,window[j])),delta=products.reduce((a,b)=>a^b,0);
  const before={sigma:[...sigma],B:[...B],gamma,L,iteration:r+1,window,delta};
  trace.push({...before,phase:'discrepancy',products});
  const scale=sigma.map(v=>f.mul(gamma,v)),correction=B.map(v=>f.mul(delta,v));
  sigma=scale.map((v,j)=>v^(j?correction[j-1]:0));
  // m_r has already incremented in ERROR_POLY_DISCREPENCY in the supplied RTL.
  if(delta!==0&&2*L<=r+1){B=[...before.sigma];L=r+1-L;gamma=delta;}
  else B=[0,...B.slice(0,4)];
  trace.push({...before,phase:'update',scale,correction,nextSigma:[...sigma],nextB:[...B],nextGamma:gamma,nextL:L});
 }
 return {trace,sigma,L};
}
export function traceWord(received){
 const n=received.length,{t}=codes[n],f=field(n),frames=[];
 let S=Array(2*t).fill(0);
 const stream=[0,...received.slice().reverse()];
 for(let q=0;q<stream.length;q+=8){
  const byte=stream.slice(q,q+8),before=[...S];
  S=S.map((v,j)=>byte.reduce((s,b)=>f.mul(s,f.pow(j+1))^b,v));
  frames.push({phase:'syndrome',byte,positions:byte.map((_,i)=>n-q-i),before,S:[...S],batch:q/8});
 }
 const synEnd=frames.length,clean=S.every(v=>v===0),bm=clean?{trace:[],sigma:[1,0,0,0,0],L:0}:inversionless(S,n);
 frames.push(...bm.trace.map(frame=>({...frame,S})));
 const bmEnd=frames.length,roots=[];let coefficients=[...bm.sigma];
 if(!clean)for(let q=0;q<n;q+=8){
  const lanes=Array.from({length:8},(_,lane)=>{
   const step=lane+1,position=n-1-q-lane,terms=coefficients.map((v,k)=>f.mul(v,f.pow(step*k))),value=terms.reduce((a,b)=>a^b,0);
   return {position,terms,value,root:position>=0&&value===0};
  });
  const next=coefficients.map((v,k)=>f.mul(v,f.pow(8*k)));
  frames.push({phase:'chien',S,sigma:[...coefficients],locator:bm.sigma,lanes,found:[...roots],nextSigma:next,L:bm.L,batch:q/8});
  for(const lane of lanes)if(lane.root)roots.push(lane.position);
  coefficients=next;
 }
 const corrected=[...received];roots.forEach(p=>corrected[p]^=1);
 const valid=roots.length===bm.L&&bm.L<=t&&syndromes(corrected).every(v=>v===0);
 const final={phase:'result',S,sigma:bm.sigma,roots,corrected,valid,clean,L:bm.L};frames.push(final);
 return {frames,S,bm,roots,corrected,valid,clean,synEnd,bmEnd,last:frames.length-1};
}
export function experiment(n=63,scenario='two',mode='hard',custom=null){
 const {k,t}=codes[n],original=encode(n,Array.from({length:k},(_,i)=>(i*7+3)%11<5?1:0));
 const presets={two:[3,17],soft:[3,...Array.from({length:t},(_,i)=>21+7*i)],many:Array.from({length:t+3},(_,i)=>5+11*i),clean:[]};
 const errors=custom??presets[scenario],received=[...original];errors.forEach(p=>received[p]^=1);
 const reliability=Array.from({length:n},(_,p)=>40+(p*13)%25);if(scenario==='soft')errors.forEach(p=>reliability[p]=8);reliability[3]=2;reliability[17]=4;
 const weak=[3,17];
 const candidates=Array.from({length:mode==='soft'?4:1},(_,mask)=>{
  const bits=[...received],flips=weak.filter((p,i)=>mask&(1<<i));flips.forEach(p=>bits[p]^=1);
  const trace=traceWord(bits),score=trace.valid?trace.corrected.reduce((sum,b,p)=>sum+(b!==received[p]?reliability[p]:0),0):Infinity;
  return {...trace,mask,bits,flips,score};
 });
 const best=candidates.filter(c=>c.valid).sort((a,b)=>a.score-b.score)[0]||null;
 return {n,original,received,errors,weak,reliability,candidates,best,recovered:!!best&&best.corrected.every((b,p)=>b===original[p])};
}
