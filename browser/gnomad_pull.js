// Run in the console of a gnomad.broadinstitute.org tab. gnomAD v4 missense positions observed in >=2 alleles = "tolerated" controls.
// Output hex -> data/variants/<GENE>.gn.hex
window.G = {};
for (const [g,N] of [['SCN1A',2010],['SCN2A',2010],['KCNQ2',872]]) {
  const q = `query{gene(gene_symbol:"${g}",reference_genome:GRCh38){canonical_transcript_id variants(dataset:gnomad_r4){variant_id consequence hgvsp exome{ac} genome{ac}}}}`;
  const j = await (await fetch('https://gnomad.broadinstitute.org/api',{method:'POST',headers:{'content-type':'application/json'},body:JSON.stringify({query:q})})).json();
  const bits = new Array(N+1).fill(0);
  for (const v of j.data.gene.variants){
    if (v.consequence!=='missense_variant' || !v.hgvsp) continue;
    const m = v.hgvsp.match(/p\.[A-Za-z]{3}(\d+)/); if(!m) continue;
    const ac = ((v.exome&&v.exome.ac)||0)+((v.genome&&v.genome.ac)||0);
    if (ac>=2) bits[+m[1]]=1;
  }
  let h=''; for(let i=0;i<=N;i+=4) h+=parseInt(bits.slice(i,i+4).join('').padEnd(4,'0'),2).toString(16);
  window.G[g] = {tx:j.data.gene.canonical_transcript_id, hex:h};
}
