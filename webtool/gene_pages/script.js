// Get the page ID from the URL 
function getQueryParam(param) {
    const urlParams = new URLSearchParams(window.location.search);
    return urlParams.get(param);
}

// Function to generate the table and append it to the page
function generateTable(pageData) {
    const container = document.getElementById('tables-container');
    
    // Check if container exists
    if (!container) {
        console.error("Table container not found!");
        return;
    }

    // Create a div for the page
    const pageDiv = document.createElement('div');
    pageDiv.classList.add('page');

    // Create table element
    const table = document.createElement('table');
    table.classList.add('table');

    // Table header
    const thead = document.createElement('thead');
    const headerRow = document.createElement('tr');
    headerRow.innerHTML = `
        <th rowspan="2" >Strata</th>
        <th rowspan="2">Detected in</th>
        <th colspan="3" class="divider-left">Best Hit in Stratum</th>
        <th colspan="6" class="divider-left">Stratum Example Species</th>
    `;
    thead.appendChild(headerRow);

    // Add sub-headers for the "Example" column
    const subHeaderRow = document.createElement('tr');
    subHeaderRow.innerHTML = `
        <th class="divider-left">Organism</th> 
        <th>E score</th>
        <th>Gene</th>
        <th class="divider-left">Organism</th> 
        <th>E score</th>
        <th>Gene</th>
        <th>Name</th>
        <th>GO terms</th>
        <th>Subcellular Loc</th>
    `;
    thead.appendChild(subHeaderRow);
    table.appendChild(thead);

    // Table body
    const tbody = document.createElement('tbody');
    const rowNames = document.querySelectorAll('#row-names li');

    pageData.rows.forEach(row => {
        const rowElement = document.createElement('tr');

        // Get the row name from the hidden list based on row number
        const rowNameCell = document.createElement('td');
        rowNameCell.textContent = rowNames[row.row - 1].textContent;  // row number - 1 to get the correct item in the list
        rowElement.appendChild(rowNameCell);

        // "Detected in" cell
        const detectedInCell = document.createElement('td');
        detectedInCell.textContent = row['detect.frac'];
        rowElement.appendChild(detectedInCell);
        

        // Best hit cells (sub-columns): best organism, E score, and gene name
        const hitsData = row.hits;
        
        const BestOrgCell = document.createElement('td');
        BestOrgCell.textContent = hitsData["BestOrg"];
        BestOrgCell.classList.add('divider-left'); // make divider
        rowElement.appendChild(BestOrgCell);
        
        const BestECell = document.createElement('td');
        BestECell.textContent = hitsData["BestE"];
        rowElement.appendChild(BestECell);

        //Make best hit into link
        const BestHitCell = document.createElement('td');

        const BestOrg = hitsData["BestOrg"]?.[0];
        const BestHit = hitsData["BestHit"]?.[0];
        const rowNum = row.row?.[0]; // get the actual row number
        const CAASOrgs = ["A632", "Huangzaosi", "Jing724","Jing92", "Dan340", "PH207", "S37", "Xu178", "Ye478", "Zheng58"]; // CAAS genomes handled separately because they don't have browsers

        let BestHitLink;
        if (CAASOrgs.includes(BestOrg)) { // CAAS genomes don't have browsers
            BestHitLink = `https://www.maizegdb.org/genome/assembly/Zm-${BestOrg}-REFERENCE-CAAS_FIL-1.0`;
        } else if (BestOrg === "Chang7 2") {
        	BestHitLink = `  https://www.maizegdb.org/genome/assembly/Zm-Chang-7_2-REFERENCE-CAAS_FIL-1.0`;
        } else if (BestOrg==="Hp301") {
        	/*BestHitLink = `https://www.maizegdb.org/gene_center/gene/${BestHit}`;*/
            BestHitLink = `https://jbrowse.maizegdb.org/?data=HP301&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`
        } else if (BestOrg==="MS71") {
        	/*BestHitLink = `https://www.maizegdb.org/gene_center/gene/${BestHit}`;*/
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Ms71&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`
        } else if (rowNum===14) {
        	/*BestHitLink = `https://www.maizegdb.org/gene_center/gene/${BestHit}`;*/
            BestHitLink = `https://jbrowse.maizegdb.org/?data=${BestOrg}&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`
        } else if (BestOrg === "Zea mays parviglumis TIL11") {
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Zv-TIL11&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (BestOrg === "Zea mays parviglumis TIL01") {
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Zv-TIL01&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (BestOrg === "Zea mays mexicana TIL25") {
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Zx-TIL25&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (BestOrg === "Zea mays mexicana TIL18") {
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Zx-TIL18&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (BestOrg === "Zea mays huehuetenagensis") {
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Zh-RIMHU001&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (BestOrg === "Zea diploperennis Gigi") {
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Zd-Gigi&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (BestOrg === "Zea diploperennis Momo") {
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Zd-Momo&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (BestOrg === "Zea nicaraguensis") {
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Zn-PI615697&loc=${BestHit}&tracks=gene_models_official%2CAv-Kellogg1287_8_gms%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (BestOrg === "Tripsacum dactyloides FL") {
            BestHitLink = `https://maizegdb.org/genome/assembly/Td-FL_9056069_6-REFERENCE-PanAnd-2.0`;
        } else if (BestOrg === "Tripsacum dactyloides KS") {
            BestHitLink = `https://maizegdb.org/genome/assembly/Td-KS_B6_1-REFERENCE-PanAnd-2.0`;
        //} else if (BestOrg === "Tripsacum dactyloides McKain") {
        //    BestHitLink = `https://maizegdb.org/genome/assembly/Td-McKain334_5-DRAFT-PanAnd-1.0`;
        } else if (BestOrg === "Anatherum virginicum") {
            BestHitLink = `https://jbrowse.maizegdb.org/?data=Av-Kellogg1287_8&loc=${BestHit}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (BestOrg === "Ginkgo biloba") {
            BestHitLink = `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies/SGTW-Ginkgo_biloba`
        } else {
            BestHitLink = `https://www.uniprot.org/uniprotkb/${BestHit}`;
        }

    
        // Create the hyperlink element
        const bestHitLinkElement  = document.createElement('a');
        bestHitLinkElement .href = BestHitLink;  // Set the href
        bestHitLinkElement .textContent = BestHit;  // Set the link text to the UniProt ID
        bestHitLinkElement .target = "_blank";  // Open in a new tab

        // Append the link to the cell
        BestHitCell.appendChild(bestHitLinkElement );
        rowElement.appendChild(BestHitCell);

        // Example hit cells (sub-columns): example organism, E score, gene name, GO terms, and subcellular location
        const ExOrgCell = document.createElement('td');
        ExOrgCell.textContent = hitsData["ExOrg"];
        ExOrgCell.classList.add('divider-left'); // make divider
        rowElement.appendChild(ExOrgCell);
        
       const ExECell = document.createElement('td');
        ExECell.textContent = hitsData["ExE"];
        rowElement.appendChild(ExECell);
        
        // Creating hyperlink for example hit
        const ExOrg = hitsData["ExOrg"]?.[0];
        const ExHit = hitsData["ExHit"]?.[0];

        let ExHitLink;
        if (ExOrg === "B73") {
        	ExHitLink = `https://jbrowse.maizegdb.org/?data=B73&loc=${encodeURIComponent(ExHit)}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (ExOrg === "Zea mays parviglumis TIL11") {
            ExHitLink = `https://jbrowse.maizegdb.org/?data=Zv-TIL11&loc=${encodeURIComponent(ExHit)}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else if (ExOrg === "Zea diploperennis Gigi") {
            ExHitLink = `https://jbrowse.maizegdb.org/?data=Zd-Gigi&loc=${encodeURIComponent(ExHit)}&tracks=gene_models_official%2Cuniprot_Escherichia_coli_090%2Cuniprot_Caenorhabditis_elegans_090%2Cuniprot_Chlamydomonas_reinhardtii_090%2Cuniprot_Chara_braunii_090%2Cuniprot_Ceratopteris_richardii_090%2Cuniprot_Arabidopsis_thaliana_090%2Cuniprot_Musa_acuminata_090%2Cuniprot_Ananas_comosus_090%2Cuniprot_Oryza_sativa_090%2Cuniprot_Panicum_virgatum_090%2Cuniprot_Sorghum_bicolor_090`;
        } else {
            ExHitLink = `https://www.uniprot.org/uniprotkb/${ExHit}`;
        }
        
        const ExHitEntryCell = document.createElement('td');

        // Create the hyperlink element
        const linkElement = document.createElement('a');
        linkElement.href = ExHitLink;  // Set the href
        linkElement.textContent = ExHit;  // Set the link text to the UniProt ID
        linkElement.target = "_blank";  // Open in a new tab

        // Append the link to the cell
        ExHitEntryCell.appendChild(linkElement);
        rowElement.appendChild(ExHitEntryCell);

        // Add protein name, GO terms, and subcellular location cells
        const proteinNameCell = document.createElement('td');
        proteinNameCell.textContent = hitsData["Protein Name"];
        rowElement.appendChild(proteinNameCell);

        const goTermsCell = document.createElement('td');
        goTermsCell.textContent = hitsData["GO terms"];
        rowElement.appendChild(goTermsCell);

        const subcellularLocCell = document.createElement('td');
        subcellularLocCell.textContent = hitsData["Subcellular Loc"];
        rowElement.appendChild(subcellularLocCell);
 
        tbody.appendChild(rowElement);
    });

    // Append the tbody to the table
    table.appendChild(tbody);

    // Append the table to the page div
    pageDiv.appendChild(table);

    // Append the page div to the container
    container.appendChild(pageDiv);
}

// When the page loads, get the page ID from the URL and generate the corresponding table. 
window.onload = function () {
    const pageId = getQueryParam('page') || window.injectedGeneID;
    console.log("Using pageId:", pageId);

    // Show spinner now
    document.getElementById("loadingOverlay").style.display = "flex";

    // load the JSON data
    fetch('cleaned_uniprot_table_loc.json')
    // fetch('https://download.maizegdb.org/data_templates/phylostrata/cleaned_uniprot_table_loc.json')
    .then(response => {
            if (!response.ok) {
                throw new Error('Failed to load JSON data');
            }
            return response.json();
        })
        .then(jsonData => {
            const pageData = jsonData.pages.find(page => page.id === pageId);
            console.log(pageData);
            if (pageData) {
                document.title = `Page: ${pageData.id}`;
                document.getElementById('page-title').textContent = pageData.id;
                generateTable(pageData);
            } else {
                console.warn('Page ID not found in JSON.');
            }
        })
        .catch(error => {
            console.error('Error loading JSON:', error);
        })
        .finally(() => {
            // Hide spinner after fetch completes (success or fail)
            document.getElementById("loadingOverlay").style.display = "none";
        });
};


