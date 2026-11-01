// Positive control 2: Nav1.7 + lamotrigine (PDB 8THH, ligand code IYJ, two sites). Vecura fpocket job b25c432e-413a-4f1d-b87d-678c9dc19c51 (chain A only).
// Ligand-contact residues (any atom within 4.5 A of the ligand, chain A), from files.rcsb.org/download/8THH.pdb:
const site1 = [391,1403,1404,1405,1406,1445,1449,1695,1696,1697,1744,1748,1752]; // HETATM A2003
const site2 = [402,403,964,967,968,1457,1461,1760];                              // HETATM A2004
// Pocket lining = protein heavy atoms within 4 A of the pocket's alpha-spheres in fpocket_pockets.pdb.
// Result (see README): site 1 -> pocket 2 (druggability 0.997) lines 11/13 contacts;
//                      site 2 -> pocket 7 (druggability 0.078) lines 7/8 contacts (found, but not rated druggable).
