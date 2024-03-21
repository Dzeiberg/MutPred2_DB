select p.ensembl_protein_id,LENGTH(p.sequence) * 19 as n_possible from proteins as p INNER JOIN (SELECT COUNT(*) as count,ensembl_protein_id from variant GROUP BY ensembl_protein_id HAVING COUNT(*) > 1) as v ON p.ensembl_protein_id = v.ensembl_protein_id AND LENGTH(p.sequence) * 19 =v.count;

