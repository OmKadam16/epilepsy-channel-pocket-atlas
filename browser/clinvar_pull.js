// Run in the console of a www.ncbi.nlm.nih.gov tab (same-origin E-utilities). Pulls ClinVar missense variants by germline class.
// Results land in window.__cv[gene] as "P:Arg553Trp"-style strings (P, LP, B, LB).
window.__cv = window.__cv || {};
const sleep = ms => new Promise(r => setTimeout(r, ms));
async function jf(u){ for (let k=0;k<6;k++){ await sleep(450+k*600); try{ return JSON.parse(await (await fetch(u)).text()); }catch(e){} } return null; }
const cl = {P:'clinsig pathogenic', LP:'clinsig likely pathogenic', B:'clinsig benign', LB:'clinsig likely benign'};
for (const g of ['SCN1A','SCN2A','KCNQ2']) {
  const out = [];
  for (const [k,c] of Object.entries(cl)) {
    const term = `${g}[gene] AND "missense variant"[molecular consequence] AND "${c}"[Properties]`;
    const s = await jf(`/entrez/eutils/esearch.fcgi?db=clinvar&retmode=json&retmax=5000&term=${encodeURIComponent(term)}`);
    const ids = s.esearchresult.idlist;
    for (let i=0;i<ids.length;i+=150){
      const d = await jf(`/entrez/eutils/esummary.fcgi?db=clinvar&retmode=json&id=${ids.slice(i,i+150).join(',')}`);
      if(!d) continue;
      for (const id of d.result.uids){
        const t = d.result[id].title, m = t.match(/\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)/);
        if (m && t.includes('('+g+')')) out.push(k+':'+m[1]+m[2]+m[3]);
      }
    }
  }
  window.__cv[g] = out;
}
// Encode unique positions as a hex bitset (bit i = residue i) -> data/variants/<GENE>.path.hex / .ben.hex
function enc(g, classes, N){ const bits=new Array(N+1).fill(0);
  for (const v of window.__cv[g]){ const [k,p]=v.split(':'); if(!classes.includes(k)) continue; bits[+p.match(/(\d+)/)[1]]=1; }
  let h=''; for(let i=0;i<=N;i+=4) h+=parseInt(bits.slice(i,i+4).join('').padEnd(4,'0'),2).toString(16); return h; }
// enc('SCN1A',['P','LP'],2010) -> pathogenic ; enc('SCN1A',['B','LB'],2010) -> benign
