// js to make blast results page, reachable from "Detected in" column of gene page
function getQueryParam(param) {
    return new URLSearchParams(window.location.search).get(param);
}

function generateBlastTable(hits, geneId, rowNum) {
    const container = document.getElementById('tables-container');

    const table = document.createElement('table');
    table.classList.add('table');

    // Header
    const thead = document.createElement('thead');
    thead.innerHTML = `
        <tr>
            <th>Organism</th>
            <th>E value</th>
            <th>Gene</th>
            <th>Score</th>
            <th>Query Start</th>
            <th>Query End</th>
            <th>Subject Start</th>
            <th>Subject End</th>
        </tr>`;
    table.appendChild(thead);

    // Body
    const tbody = document.createElement('tbody');
    hits.forEach(hit => {
    const tr = document.createElement('tr');

    const BestOrg = hit.organism;
    const BestHit = hit.sseqid;
    const rowNum = hit.phylostratum;

    const CAASOrgs = ["A632", "Huangzaosi", "Jing724", "Jing92", "Dan340", "PH207", "S37", "Xu178", "Ye478", "Zheng58"];

    let hitLink;
    if (CAASOrgs.includes(BestOrg)) {
        hitLink = `https://www.maizegdb.org/genome/assembly/Zm-${BestOrg}-REFERENCE-CAAS_FIL-1.0`;
    } else if (BestOrg === "Chang7 2") {
        hitLink = `https://www.maizegdb.org/genome/assembly/Zm-Chang-7_2-REFERENCE-CAAS_FIL-1.0`;
    } else if (BestOrg === "Hp301") {
        hitLink = `https://jbrowse.maizegdb.org/?data=HP301&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
    } else if (BestOrg === "MS71") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Ms71&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
    } else if (rowNum === 14) {
        hitLink = `https://jbrowse.maizegdb.org/?data=${BestOrg}&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
    } else if (BestOrg === "Zea mays parviglumis TIL11") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Zv-TIL11&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
    } else if (BestOrg === "Zea mays parviglumis TIL01") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Zv-TIL01&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
    } else if (BestOrg === "Zea mays mexicana TIL25") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Zx-TIL25&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_₀`;
    } else if (BestOrg === "Zea mays mexicana TIL18") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Zx-TIL18&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
    } else if (BestOrg === "Zea mays huehuetenagensis") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Zh-RIMHU001&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_₀`;
    } else if (BestOrg === "Zea diploperennis Gigi") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Zd-Gigi&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_₀`;
    } else if (BestOrg === "Zea diploperennis Momo") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Zd-Momo&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_09₀`;
    } else if (BestOrg === "Zea nicaraguensis") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Zn-PI615697&loc=${BestHit}&tracks=gene_models_official%2CAv-Kellogg1287_8_gms%2Cuniprot_Escherichia_coli_₀`;
        hitLink = `https://maizegdb.org/genome/assembly/Td-FL_9056069_6-REFERENCE-PanAnd-2.0`;
    } else if (BestOrg === "Tripsacum dactyloides KS") {
        hitLink = `https://maizegdb.org/genome/assembly/Td-KS_B6_1-REFERENCE-PanAnd-2.0`;
    } else if (BestOrg === "Anatherum virginicum") {
        hitLink = `https://jbrowse.maizegdb.org/?data=Av-Kellogg1287_8&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
    } else if (BestOrg === "Ginkgo biloba") {
        hitLink = `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies/SGTW-Ginkgo_biloba`;
    } else {
        hitLink = `https://www.uniprot.org/uniprotkb/${BestHit}`;
    }

    tr.innerHTML = `
        <td>${hit.organism}</td>
        <td>${hit.evalue}</td>
        <td><a href="${hitLink}" target="_blank">${hit.sseqid}</a></td>
        <td>${hit.score}</td>
        <td>${hit.qstart}</td>
        <td>${hit.qend}</td>
        <td>${hit.sstart}</td>
        <td>${hit.send}</td>
    `;
    tbody.appendChild(tr);
});

    table.appendChild(tbody);
    container.appendChild(table);
}

window.onload = function () {
    const geneId = getQueryParam('page');
    const rowNum = parseInt(getQueryParam('row'));

    // Set back link
    document.getElementById('back-link').href = `index.html?page=${encodeURIComponent(geneId)}`;
    document.getElementById('page-title').textContent =
        `BLAST Results — ${geneId} (Search Space: Phylostratum ${rowNum})`;
    document.title = `BLAST: ${geneId} PS${rowNum}`;

    document.getElementById('loadingOverlay').style.display = 'flex';


    // NOTE: one example json of BLAST results is provided. They should be generated for all genes and added to the gene_pages/blast_results_by_gene/ directory as {geneId}_blast_results.json
    fetch(`blast_results_by_gene/${encodeURIComponent(geneId)}_blast_results.json`)
        .then(r => {
            if (!r.ok) throw new Error('Failed to load BLAST data');
            return r.json();
        })
        .then(data => {
            const hits = data.blast_results.filter(
                d => d.id === geneId && d.phylostratum === rowNum
            );
            if (hits.length > 0) {
                generateBlastTable(hits, geneId, rowNum);
            } else {
                document.getElementById('tables-container').textContent =
                    'No BLAST results found for this gene and stratum.';
            }
        })
        .catch(err => console.error('Error:', err))
        .finally(() => {
            document.getElementById('loadingOverlay').style.display = 'none';
        });
};