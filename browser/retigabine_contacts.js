// Positive control: residues within 4.5 A of retigabine (ligand code FBX) in KCNQ2 7CR2. Run on any page that can fetch files.rcsb.org.
const t = await (await fetch('https://files.rcsb.org/download/7CR2.pdb')).text(), L = t.split('\n');
const lig = L.filter(l=>l.startsWith('HETATM')&&l.substr(17,3)==='FBX').map(l=>({c:l[21],x:+l.substr(30,8),y:+l.substr(38,8),z:+l.substr(46,8)}));
const at = L.filter(l=>l.startsWith('ATOM')&&l[77]!=='H').map(l=>({c:l[21],n:+l.substr(22,4),x:+l.substr(30,8),y:+l.substr(38,8),z:+l.substr(46,8)}));
// ... residues with any atom within 4.5 A of any ligand atom (chain A): 104,236,239,240,243,272,299,300,303,304,305,308,312
