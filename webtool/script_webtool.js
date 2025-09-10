let downloadCanceled = false;

let originalGeneOrder = []; // Store the original input order
let geneDataCache = null; // Cache the gene data to avoid repeated fetches
let currentSortDirection = 'asc'; // Track current sort direction

// check if genes are valid
let validGeneSet = new Set();
let jsonReady = false;

fetch('image_code/gene_fill.json')
  .then(response => response.json())
  .then(data => {
    data.pages.forEach(entry => {
      validGeneSet.add(entry.id);
    });
    jsonReady = true;
    console.log("Valid gene set loaded:", validGeneSet);
  })
  .catch(error => {
    console.error("Error loading gene JSON:", error);
  });

let cleanGeneIds = []; // Store cleaned gene IDs after validation

function handleSearch() {
  const input = document.getElementById("geneInput").value.trim();
  const geneIds = input.split(/\s*\n\s*/).filter(id => id !== "");
  
  // Limit to 100 gene IDs
  if (geneIds.length > 100) {
    alert("You have entered more than 100 genes.\n\nPlease submit fewer genes at one time, or use the Downloads page to download information for the whole genome.");
    return;
  }

  // And clear the gene data cache when starting a new search:
  geneDataCache = null; // Reset cache for new search

  // Validate gene IDs and find invalid and duplicate IDs
  const seen = new Set();
  cleanGeneIds = [];
  const missingGenes = [];
  const duplicateGenes = [];

  geneIds.forEach(id => {
    if (!validGeneSet.has(id)) {
      missingGenes.push(id);
      return; // skip invalid
    }
    
    if (seen.has(id)) {
      if (!duplicateGenes.includes(id)) {
        duplicateGenes.push(id);
      }
      return; // skip duplicate
    }
    
    seen.add(id);
    cleanGeneIds.push(id); // add only first, valid occurrence
  });

  if(cleanGeneIds.length === 0) {
    alert("No valid gene IDs found. Please check your input.");
    return;
  }

  originalGeneOrder = [...cleanGeneIds]; // Store original input order

  // Also reset the sort dropdown to "input" and direction to ascending on new searches:
  document.getElementById("sortSelect").value = "input";
  currentSortDirection = 'asc';
  updateSortDirectionButton();

  // Now use only cleanGeneIds to generate joinedIds
  const joinedIds = cleanGeneIds.join(",");
  const iframeUrl = `image_code/index.html?page=${encodeURIComponent(joinedIds)}`;

  //insert iframe for search result image
  const container = document.getElementById("imageResults");
  container.innerHTML = ""; // Clear previous

  // After you've created cleanGeneIds, use to set dimensions of image
  const baseHeight = 300;      // starting height in px
  const heightPerGene = 50;    // additional height per gene

  const iframeHeight = baseHeight + heightPerGene * cleanGeneIds.length;

  const iframe = document.createElement("iframe");
  iframe.src = iframeUrl;
  iframe.width = "1250";         
  iframe.height = iframeHeight;  // dynamic height
  iframe.id = "dynamicImageIframe";
  iframe.style.border = "none";
  iframe.style.display = "block";
  iframe.style.margin = "0 auto";

  container.appendChild(iframe);

  const warningBox = document.getElementById("missingGenesWarning");

  warningBox.style.display = "none";
  warningBox.innerText = "";

  const resultsBar = document.getElementById("searchResultsBar");
  resultsBar.style.display = geneIds.length > 0 ? "flex" : "none";

  document.getElementById("selectAllCheckbox").checked = false; // Reset checkboxes on new search

  updateWarningBox(missingGenes, duplicateGenes); // Show warnings after search loop ends

  function updateWarningBox(missing, duplicates) {
    warningBox.innerHTML = ""; // Clear old content

    if (missing.length > 0) {
      const warningDiv = document.createElement("div");
      warningDiv.className = "warning-message";
      warningDiv.textContent = `Warning: No results found for: ${missing.join(", ")}`;
      warningBox.appendChild(warningDiv); // Append to warningBox
    }

    if (duplicates.length > 0) {
      const noteDiv = document.createElement("div");
      noteDiv.className = "note-message";
      noteDiv.textContent = `Note: Duplicate searches ignored for: ${duplicates.join(", ")}`;
      warningBox.appendChild(noteDiv); // Append to warningBox
    }

    if (missing.length > 0 || duplicates.length > 0) {
      warningBox.style.display = "block";
    }
  }

}

// checkboxes
let selectedGenesFromIframe = [];
window.addEventListener('message', function(event) {
    if (event.origin !== window.location.origin) return;

    if (event.data.type === 'updateGeneSelection') {
        selectedGenesFromIframe = event.data.selectedGenes;
        
        // Update UI to show selection count
        document.getElementById("selectedCount").textContent = selectedGenesFromIframe.length;
        
        // Update the "Select All" checkbox state based on iframe checkboxes
        const selectAllCheckbox = document.getElementById("selectAllCheckbox");
        if (event.data.allSelected) {
            selectAllCheckbox.checked = true;
            selectAllCheckbox.indeterminate = false;
        } else if (event.data.noneSelected) {
            selectAllCheckbox.checked = false;
            selectAllCheckbox.indeterminate = false;
        } else {
            // Some but not all are selected - show indeterminate state
            selectAllCheckbox.checked = false;
            selectAllCheckbox.indeterminate = true;
        }
    }
});

// load gene data for sorting
async function loadGeneDataForSorting() {
  if (geneDataCache) {
    return geneDataCache; // Return cached data if available
  }
  
  try {
    const response = await fetch('image_code/gene_fill.json');
    const data = await response.json();
    geneDataCache = data;
    return data;
  } catch (error) {
    console.error("Error loading gene data for sorting:", error);
    return null;
  }
}

// sort genes based on the selected option
async function sortGenes(sortOption, direction = 'asc') {
  if (cleanGeneIds.length === 0) return;
  
  let sortedIds = [...cleanGeneIds]; // Create a copy
  
  switch (sortOption) {
    case 'alphabetic':
      sortedIds.sort((a, b) => {
        const result = a.localeCompare(b);
        return direction === 'asc' ? result : -result;
      });
      break;
      
    case 'phylostrata':
      const geneData = await loadGeneDataForSorting();
      if (geneData) {
        // Create a map of gene ID to fill_pct for quick lookup
        const fillPctMap = {};
        geneData.pages.forEach(page => {
          fillPctMap[page.id] = page.fill_pct;
        });
        
        // Sort by fill_pct (phylostrata), with fallback to alphabetic for missing data
        sortedIds.sort((a, b) => {
          const fillA = fillPctMap[a] || 0;
          const fillB = fillPctMap[b] || 0;
          
          if (fillA !== fillB) {
            const result = fillA - fillB;
            return direction === 'asc' ? result : -result;
          }
          const alphaResult = a.localeCompare(b); // Fallback to alphabetic
          return direction === 'asc' ? alphaResult : -alphaResult;
        });
      }
      break;
      
    case 'input':
    default:
      sortedIds = [...originalGeneOrder]; // Start with original order
      if (direction === 'desc') {
        sortedIds.reverse(); // Reverse for descending
      }
      break;
  }
  
  // Update cleanGeneIds with the new order
  cleanGeneIds = sortedIds;
  
  // Regenerate the iframe with the new order
  regenerateIframe();
}

// regenerate the iframe
function regenerateIframe() {
  const joinedIds = cleanGeneIds.join(",");
  const iframeUrl = `image_code/index.html?page=${encodeURIComponent(joinedIds)}`;

  // Regenerate iframe with new order
  const container = document.getElementById("imageResults");
  container.innerHTML = ""; // Clear previous

  const baseHeight = 300;
  const heightPerGene = 50;
  const iframeHeight = baseHeight + heightPerGene * cleanGeneIds.length;

  const iframe = document.createElement("iframe");
  iframe.src = iframeUrl;
  iframe.width = "1250";
  iframe.height = iframeHeight;
  iframe.id = "dynamicImageIframe";
  iframe.style.border = "none";
  iframe.style.display = "block";
  iframe.style.margin = "0 auto";

  container.appendChild(iframe);
  
  // Reset the select all checkbox since we have a new iframe
  document.getElementById("selectAllCheckbox").checked = false;
  document.getElementById("selectAllCheckbox").indeterminate = false;
}


// sort event handler
function handleSortChange() {
  const sortSelect = document.getElementById("sortSelect");
  const selectedSort = sortSelect.value;
  
  // Show loading state
  const container = document.getElementById("imageResults");
  const currentIframe = container.querySelector("iframe");
  if (currentIframe) {
    currentIframe.style.opacity = "0.5";
  }
  
  // Sort and regenerate
  sortGenes(selectedSort, currentSortDirection).then(() => {
    // Remove loading state
    const newIframe = container.querySelector("iframe");
    if (newIframe) {
      newIframe.style.opacity = "1";
    }
  }).catch(error => {
    console.error("Error sorting genes:", error);
    // Remove loading state even on error
    const iframe = container.querySelector("iframe");
    if (iframe) {
      iframe.style.opacity = "1";
    }
  });
}

function toggleSortDirection() {
  // Toggle direction
  currentSortDirection = currentSortDirection === 'asc' ? 'desc' : 'asc';
  
  // Update button appearance
  updateSortDirectionButton();
  
  // Re-sort with new direction
  const sortSelect = document.getElementById("sortSelect");
  const selectedSort = sortSelect.value;
  
  // Show loading state
  const container = document.getElementById("imageResults");
  const currentIframe = container.querySelector("iframe");
  if (currentIframe) {
    currentIframe.style.opacity = "0.5";
  }
  
  // Sort and regenerate
  sortGenes(selectedSort, currentSortDirection).then(() => {
    // Remove loading state
    const newIframe = container.querySelector("iframe");
    if (newIframe) {
      newIframe.style.opacity = "1";
    }
  }).catch(error => {
    console.error("Error sorting genes:", error);
    // Remove loading state even on error
    const iframe = container.querySelector("iframe");
    if (iframe) {
      iframe.style.opacity = "1";
    }
  });
}

function updateSortDirectionButton() {
  const btn = document.getElementById("sortDirectionBtn");
  if (currentSortDirection === 'asc') {
    btn.innerHTML = "↑ Asc";
    btn.title = "Click to sort descending";
  } else {
    btn.innerHTML = "↓ Desc";
    btn.title = "Click to sort ascending";
  }
}

// Move among the different sections of the webtool
function showSectionFromHash() {
  const hash = window.location.hash || "#search";
  document.getElementById("search")?.style?.setProperty("display", "none");
  document.getElementById("downloadsSection")?.style?.setProperty("display", "none");
  document.getElementById("aboutSection")?.style?.setProperty("display", "none");

  if (hash === "#about") {
    document.getElementById("aboutSection").style.display = "block";
    // loadAboutPage(); // Only when needed
  } else if (hash === "#downloads") {
    document.getElementById("downloadsSection").style.display = "block";
  } else {
    document.getElementById("search").style.display = "block";
  }
}

window.addEventListener("hashchange", showSectionFromHash);
window.addEventListener("DOMContentLoaded", showSectionFromHash);

function toggleSelectAll(checkbox) {
  // Send message to iframe to toggle all checkboxes
  const iframe = document.getElementById("dynamicImageIframe");
  if (iframe && iframe.contentWindow) {
    iframe.contentWindow.postMessage({
      type: 'toggleAllCheckboxes',
      checked: checkbox.checked
    }, '*');
  }
}

function showProgressBar() {
  document.getElementById("loadingOverlay").style.display = "flex";
  document.getElementById("cancelButton").style.display = "block";
  downloadCanceled = false;
  updateProgressBar(0);
}

function hideProgressBar() {
  document.getElementById("loadingOverlay").style.display = "none";
  document.getElementById("cancelButton").style.display = "none";
}

function cancelDownload() {
  downloadCanceled = true;
}

function updateProgressBar(percent) {
  const fill = document.getElementById("progressBarFill");
  fill.style.width = `${percent}%`;
  fill.textContent = `${percent}%`;
}

  function populateExampleGenes() {
    const exampleGenes = [
      // "Zm00001eb181160",
      // "Zm00001eb071070",
      // "Zm00001eb079890",
      // "Zm00001eb416250"
      "Zm00001eb175160",
      "Zm00001eb267570",
      "Zm00001eb281380",
      "Zm00001eb286860",
      "Zm00001eb002510",
      "Zm00001eb003760"
    ];
    document.getElementById("geneInput").value = exampleGenes.join("\n");
    handleSearch(); // Automatically trigger the search
  }

function openDownloadModal() {
  document.getElementById("allCount").textContent = cleanGeneIds.length;

  document.getElementById("downloadModal").style.display = "flex";
}

function closeDownloadModal() {
  document.getElementById("downloadModal").style.display = "none";
}

function handleDownload() {
  const scope = document.querySelector('input[name="downloadScope"]:checked').value;
  const format = document.querySelector('input[name="downloadFormat"]:checked').value;

  if (format === "images") {
    if (scope === "selected") {
      downloadSelectedGenesAsPng(scope);
    } else {
      // downloadIframeAsPng();
      downloadSelectedGenesAsPng(scope);
    }
  } else if (format === "details") {
    downloadDetails(scope);
  } else if (format === "details-text") {  // text format option
    downloadText(scope);
  }

  closeDownloadModal();
}

// fetch gene details HTML
async function fetchGeneDetailHTML(geneID) {
  // Fetch base index.html with the injected geneID script + base tag (same as before)
  const response = await fetch(`gene_pages/index.html?page=${geneID}`);
  if (!response.ok) throw new Error(`Failed to fetch page for ${geneID}`);
  
  let html = await response.text();

  const baseTag = '<base href="/gene_pages/">';
  const fallbackScript = `<script>window.injectedGeneID = "${geneID}";</script>`;

  html = html.replace(
    /<head[^>]*>/i,
    match => `${match}\n    ${baseTag}\n    ${fallbackScript}`
  );

  // create an offscreen iframe to render and wait for content 
  return new Promise((resolve, reject) => {
    const iframe = document.createElement('iframe');
    iframe.style.display = 'none';
    document.body.appendChild(iframe);

    const doc = iframe.contentWindow.document;
    doc.open();
    doc.write(html);
    doc.close();

    const maxWait = 10000;
    const pollInterval = 200;
    let waited = 0;

    const checkReady = setInterval(() => {
      const innerDoc = iframe.contentDocument || iframe.contentWindow.document;
      const isReady = innerDoc.querySelector("#tables-container .page");

      if (isReady || waited >= maxWait) {
        clearInterval(checkReady);

        if (!isReady) {
          console.warn(`Timeout reached. Page may not be fully rendered for gene ${geneID}. Proceeding with partial HTML.`);
        }

        const fullHTML = '<!DOCTYPE html>\n' + innerDoc.documentElement.outerHTML;

        document.body.removeChild(iframe);

        resolve(fullHTML);
      }

      waited += pollInterval;
    }, pollInterval);
  });
}

// functions to download text versions of detail files
function generateTextTable(pageData) {
    // Create row name list 
    const rowNames = [
        "1: Cellular organisms",
        "2: Domain: Eukaryota",
        "3: Kingdom: Viridiplantae",
        "4: Phylum: Streptophyta",
        "5: Land plants (Embryophyta)",
        "6: Class: Magnoliopsida",
        "7: Monocots (Liliopsida)",
        "8: Order: Poales",
        "9: Family: Poaceae",
        "10: Subfamily: Panicoideae",
        "11: Tribe: Andropogoneae",
        "12: Genus: Zea",
        "13: Species: Zea mays",
        "14: Subspecies: Zea mays mays"
    ];

    let text = "";
    text += `Strata\tdetect.frac\tBestOrg\tBestE\tBestHit\tExOrg\tExE\tExHit\tName\tGO terms\tSubcellular Loc\n`;

    pageData.rows.forEach(row => {
        const rowName = rowNames[row.row - 1] || "";
        const detectFrac = row['detect.frac'];
        const hits = row.hits;

        const BestOrg = hits["BestOrg"] || "";
        const BestE = hits["BestE"] || "";
        const BestHit = hits["BestHit"] || "";
        const ExOrg = hits["ExOrg"] || "";
        const ExE = hits["ExE"] || "";
        const ExHit = hits["ExHit"] || "";
        const proteinName = hits["Protein Name"] || "";
        const goTerms = hits["GO terms"] || "";
        const subLoc = hits["Subcellular Loc"] || "";

        text += `${rowName}\t${detectFrac}\t${BestOrg}\t${BestE}\t${BestHit}\t${ExOrg}\t${ExE}\t${ExHit}\t${proteinName}\t${goTerms}\t${subLoc}\n`;
    });

    return text;
}

function downloadTextFile(filename, textContent) {
    const blob = new Blob([textContent], { type: "text/plain" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

// Unified download text function
async function downloadText(scope) {
  // Determine which genes to process based on scope
  let genesToProcess;
  
  if (scope === "selected") {
    if (selectedGenesFromIframe.length === 0) {
      alert("Please select at least one gene to download.");
      return;
    }
    genesToProcess = selectedGenesFromIframe.map(gene => gene.geneId);
  } else {
    // scope === "all"
    genesToProcess = cleanGeneIds;
    if (genesToProcess.length === 0) {
      alert("No genes available to download.");
      return;
    }
  }

  // Handle single file case (no zip needed)
  if (genesToProcess.length === 1) {
    const geneID = genesToProcess[0];
    
    try {
      const response = await fetch(`gene_pages/cleaned_uniprot_table_loc.json`);
      // const response = await fetch('https://download.maizegdb.org/data_templates/phylostrata/cleaned_uniprot_table_loc.json')

      if (!response.ok) throw new Error("Failed to load JSON");
      const jsonData = await response.json();
      
      const pageData = jsonData.pages.find(p => p.id === geneID);
      if (!pageData) {
        alert(`Gene ${geneID} not found in JSON.`);
        return;
      }
      
      const textContent = generateTextTable(pageData);
      downloadTextFile(`${geneID}.txt`, textContent);
    } catch (err) {
      console.error(err);
    }
    return;
  }

  // Handle multiple files (zip creation)
  const total = genesToProcess.length;
  let completed = 0;

  showProgressBar();
  const zip = new JSZip();

  try {
    const response = await fetch(`gene_pages/cleaned_uniprot_table_loc.json`);
    // const response = await fetch('https://download.maizegdb.org/data_templates/phylostrata/cleaned_uniprot_table_loc.json')

    if (!response.ok) throw new Error("Failed to load JSON");
    const jsonData = await response.json();

    for (const geneID of genesToProcess) {
      if (downloadCanceled) {
        console.warn("Download canceled by user.");
        hideProgressBar();
        return;
      }

      const pageData = jsonData.pages.find(p => p.id === geneID);

      if (pageData) {
        const textContent = generateTextTable(pageData);
        zip.file(`${geneID}.txt`, textContent);
      } else {
        console.warn(`Gene ${geneID} not found in JSON.`);
      }

      completed++;
      const percent = Math.round((completed / total) * 100);
      updateProgressBar(percent);
    }

    const content = await zip.generateAsync({ type: "blob" });
    saveAs(content, "Phylostrata_details_text.zip");
    updateProgressBar(100);
    setTimeout(hideProgressBar, 800);

  } catch (error) {
    console.error(error);
    hideProgressBar();
  }
}

// Unified download details function
async function downloadDetails(scope) {
  // Determine which genes to process based on scope
  let genesToProcess;
  
  if (scope === "selected") {
    if (selectedGenesFromIframe.length === 0) {
      alert("Please select at least one gene detail to download.");
      return;
    }
    genesToProcess = selectedGenesFromIframe.map(gene => gene.geneId);
  } else {
    // scope === "all"
    genesToProcess = cleanGeneIds;
    if (genesToProcess.length === 0) {
      alert("No gene details available to download.");
      return;
    }
  }

  // Handle single file case (no zip needed)
  if (genesToProcess.length === 1) {
    const geneID = genesToProcess[0];
    
    try {
      showProgressBar();
      const htmlContent = await fetchGeneDetailHTML(geneID);
      const blob = new Blob([htmlContent], { type: "text/html" });
      saveAs(blob, `${geneID}.html`);
      updateProgressBar(100);
      setTimeout(hideProgressBar, 800);
    } catch (e) {
      console.error(`Failed to download detail for ${geneID}`, e);
      hideProgressBar();
    }
    return;
  }

  // Handle multiple files (zip creation)
  const total = genesToProcess.length;
  let completed = 0;

  showProgressBar();
  const zip = new JSZip();

  for (const geneID of genesToProcess) {
    if (downloadCanceled) {
      console.warn("Download canceled by user.");
      hideProgressBar();
      return;
    }

    try {
      const htmlContent = await fetchGeneDetailHTML(geneID);
      zip.file(`${geneID}.html`, htmlContent);
    } catch (e) {
      console.error(`Failed to fetch gene detail for ${geneID}:`, e);
    }

    completed++;
    const percent = Math.round((completed / total) * 100);
    updateProgressBar(percent);
  }

  zip.generateAsync({ type: "blob" }).then(content => {
    saveAs(content, "Phylostrata_details_html.zip");
    updateProgressBar(100);
    setTimeout(hideProgressBar, 800);
  });
}

// The following functions enable downloading results image as PNG
async function downloadSelectedGenesAsPng(scope) {
  let genesToProcess;

  if (scope === "selected") {
    if (selectedGenesFromIframe.length === 0) {
      alert("Please select at least one gene detail to download.");
      return;
    }
    genesToProcess = selectedGenesFromIframe.map(gene => gene.geneId);
  } else {
    // scope === "all"
    genesToProcess = cleanGeneIds;
    if (genesToProcess.length === 0) {
      alert("No gene details available to download.");
      return;
    }
  }

  try {
    showProgressBar();
    updateProgressBar(10);
    
    // Create a temporary iframe for the selected genes
    const tempIframe = await createTempIframeForSelectedGenes(genesToProcess);
    updateProgressBar(30);
    
    // Wait for the iframe to fully load
    await waitForIframeLoad(tempIframe);
    updateProgressBar(50);
    
    // Convert images and download
    const iframeDoc = tempIframe.contentDocument || tempIframe.contentWindow.document;
    
    
    if (!iframeDoc) {
      console.error('ERROR: Cannot access iframe document');
      throw new Error("Cannot access temporary iframe content");
    }
    

    const checkboxes = iframeDoc.querySelectorAll('.img-checkbox');
    checkboxes.forEach((checkbox, index) => {
      const styles = iframeDoc.defaultView.getComputedStyle(checkbox);

    });
    
    // Convert all images to data URLs first
    await convertAllImagesToDataUrls(iframeDoc, tempIframe);
    updateProgressBar(70);
    
    // Add extra wait time to ensure everything is rendered
    await new Promise(resolve => setTimeout(resolve, 2000));
    updateProgressBar(80);
    
    // Capture the entire body
    const bodyElement = iframeDoc.body;
    
    console.log('Starting html2canvas capture of selected genes...');

    const canvas = await html2canvas(bodyElement, {
      useCORS: true,
      allowTaint: true,
      scale: 3,
      logging: true, // Enable logging for debugging
      height: bodyElement.scrollHeight,
      width: bodyElement.scrollWidth,
      scrollX: 0,
      scrollY: 0,
      imageTimeout: 0,
      onclone: (clonedDoc) => {        
        // Check if checkboxes exist in cloned document
        const clonedCheckboxes = clonedDoc.querySelectorAll('.img-checkbox');
        // console.log('Checkboxes in cloned doc:', clonedCheckboxes.length);
        
        // Apply multiple methods to hide checkboxes
        const style = clonedDoc.createElement('style');
        style.textContent = `
          * { overflow: visible !important; }
          body { min-height: 800px !important; padding-bottom: 100px !important; }
          .timeline-tickmark::before { white-space: pre-wrap !important; }
          .img-checkbox { 
            display: none !important; 
            visibility: hidden !important;
            opacity: 0 !important;
            width: 0 !important;
            height: 0 !important;
          }
          input[type="checkbox"] {
            display: none !important;
            visibility: hidden !important;
          }
        `;
        clonedDoc.head.appendChild(style);
        
        // Also directly hide checkboxes via DOM manipulation
        clonedCheckboxes.forEach(checkbox => {
          checkbox.style.display = 'none';
          checkbox.style.visibility = 'hidden';
          checkbox.style.opacity = '0';
          checkbox.remove(); // Remove entirely
        });
        
        // console.log('Applied checkbox hiding styles and removed elements');
      }
    });
    
    updateProgressBar(90);
    
    // Create and download the image
    canvas.toBlob(function(blob) {
      if (blob) {
        const timestamp = new Date().toISOString().slice(0, 19).replace(/:/g, '-');
        const filename = `Phylostrata_search_result.png`;
        saveAs(blob, filename);
        updateProgressBar(100);
        setTimeout(() => {
          hideProgressBar();
          // Clean up temporary iframe
          if (document.body.contains(tempIframe)) {
            document.body.removeChild(tempIframe);
          }
        }, 800);
      } else {
        throw new Error('Failed to create image file');
      }
    }, 'image/png', 0.8);
    
  } catch (error) {
    console.error('Download failed:', error);
    alert('Download failed: ' + error.message);
    hideProgressBar();
    // Clean up any temporary iframe that might exist
    const tempIframes = document.querySelectorAll('.temp-download-iframe');
    tempIframes.forEach(iframe => {
      if (document.body.contains(iframe)) {
        document.body.removeChild(iframe);
      }
    });
  }
}

function createTempIframeForSelectedGenes(selectedGenes) {
  return new Promise((resolve) => {
    const joinedIds = selectedGenes.join(",");
    const iframeUrl = `image_code/index.html?page=${encodeURIComponent(joinedIds)}&hideCheckboxes=true`;

    const baseHeight = 300;
    const heightPerGene = 50;
    const iframeHeight = baseHeight + heightPerGene * selectedGenes.length;

    const tempIframe = document.createElement("iframe");
    tempIframe.src = iframeUrl;
    tempIframe.width = "1250";
    tempIframe.height = iframeHeight;
    tempIframe.className = "temp-download-iframe";
    tempIframe.style.position = "absolute";
    tempIframe.style.top = "-9999px";
    tempIframe.style.left = "-9999px";
    tempIframe.style.border = "none";
    tempIframe.style.display = "block";

    tempIframe.onload = () => {
      console.log('Temp iframe loaded.'); // DEBUG
      resolve(tempIframe);
    };

    document.body.appendChild(tempIframe);
  });
}

function waitForIframeLoad(iframe) {
  return new Promise((resolve) => {
    if (iframe.contentDocument && iframe.contentDocument.readyState === 'complete') {
      // Already loaded
      setTimeout(resolve, 500); // Small delay to ensure everything is rendered
    } else {
      // Wait for load event
      iframe.onload = () => {
        setTimeout(resolve, 500); // Small delay to ensure everything is rendered
      };
    }
  });
}

async function convertAllImagesToDataUrls(doc, iframe) {
  const images = Array.from(doc.querySelectorAll('img'));
  
  for (let img of images) {
    if (!img.src.startsWith('data:')) {
      try {
        const dataUrl = await simpleImageToDataUrl(img.src, iframe);
        if (dataUrl) {
          img.src = dataUrl;
        }
      } catch (error) {
        console.warn('Failed to convert image:', error);
      }
    }
  }
}

function simpleImageToDataUrl(src, iframe) {
  return new Promise((resolve) => {
    const img = new Image();
    img.crossOrigin = 'anonymous';
    
    img.onload = () => {
      try {
        const canvas = document.createElement('canvas');
        canvas.width = img.naturalWidth;
        canvas.height = img.naturalHeight;
        canvas.getContext('2d').drawImage(img, 0, 0);
        resolve(canvas.toDataURL('image/png'));
      } catch (e) {
        resolve(null);
      }
    };
    
    img.onerror = () => resolve(null);
    
    // Build URL
    if (src.startsWith('http')) {
      img.src = src;
    } else if (src.startsWith('/')) {
      img.src = window.location.origin + src;
    } else {
      const base = new URL(iframe.src).href.replace(/\/[^\/]*$/, '/');
      img.src = new URL(src, base).href;
    }
  });
}

// download all NAM genome-wide files
async function downloadZipByPattern(pattern, zipName) {
  const overlay = document.getElementById('loadingOverlayDownloads');
  overlay.style.display = 'flex'; // Show loading overlay

  const zip = new JSZip();
  const links = document.querySelectorAll('a[download]');
  const folder = zip.folder(pattern);

  try {
    for (const link of links) {
      const href = link.getAttribute('href');
      if (href.includes(pattern)) {
        const filename = href.split('/').pop();
        const response = await fetch(href);
        if (!response.ok) throw new Error(`Failed to fetch ${filename}`);
        const blob = await response.blob();
        folder.file(filename, blob);
      }
    }

    const zipBlob = await zip.generateAsync({ type: 'blob' });
    saveAs(zipBlob, zipName);
  } catch (error) {
    console.error('Error during ZIP creation:', error);
    alert('There was a problem downloading the files. Please try again.');
  } finally {
    overlay.style.display = 'none'; // Always hide overlay when done
  }
}

document.getElementById('download-details-link').addEventListener('click', (e) => {
  e.preventDefault();
  downloadZipByPattern('phylostrata_details', 'NAM_phylostrata_details.zip');
});

document.getElementById('download-results-link').addEventListener('click', (e) => {
  e.preventDefault();
  downloadZipByPattern('phylostrata_results', 'NAM_phylostrata_results.zip');
});



