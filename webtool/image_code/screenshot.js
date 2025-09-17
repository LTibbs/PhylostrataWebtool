// This code enables bulk generation of single-gene results images
// to use for MaizeGDB gene pages.
const puppeteer = require('puppeteer');
const fs = require('fs');
const path = require('path');

(async () => {
    const data = JSON.parse(fs.readFileSync('image_code/gene_fill.json', 'utf8'));
    const ids = data.pages.map(p => p.id);
    console.log(`Found ${ids.length} genes.`);

    const outputDir = 'output_images';
    if (!fs.existsSync(outputDir)) fs.mkdirSync(outputDir);

    const browser = await puppeteer.launch();
    const page = await browser.newPage();

    await page.setViewport({
        width: 1200,
        height: 350,  // give enough height to fit content
        deviceScaleFactor: 2
    });

// Loop through each gene ID to save result image
for (const id of ids) {
    const screenshotPath = path.join(outputDir, `${id}.png`);
    console.log(`Processing ${id}...`);
    // Skip if file already exists
    if (fs.existsSync(screenshotPath)) {
        console.log(`Skipping ${id}, image already exists.`);
        continue;
    }

    const url = `http://localhost:8080/image_code/index.html?page=${id}&hideCheckboxes=true`;
    // const url = `http://localhost:8080/index.html?page=${id}`;
    await page.goto(url, { waitUntil: 'networkidle0' });

    await page.evaluate(() => {
        
 });

    // Wait to ensure everything is rendered
await new Promise(resolve => setTimeout(resolve, 1000));


    const clip = await page.evaluate(() => {
        function getBox(selector) {
            const el = document.querySelector(selector);
            if (!el) return null;
            const rect = el.getBoundingClientRect();
            return {
                top: rect.top,
                left: rect.left,
                right: rect.right,
                bottom: rect.bottom
            };
        }

        const title = getBox('#page-title');
        const thermometer = getBox('.thermometer-container');
        const timeline = getBox('.timeline');

        if (!title || !thermometer || !timeline) return null;

        const top = Math.min(title.top, thermometer.top, timeline.top);
        const left = Math.min(title.left, thermometer.left, timeline.left);
        const right = Math.max(title.right, thermometer.right, timeline.right);
        const bottom = Math.max(title.bottom, thermometer.bottom, timeline.bottom + 80);  // add buffer for labels

        const padding = 20; // adjust as needed

        return {
            x: Math.max(0, left - padding),
            y: Math.max(0, top - padding),
            width: (right - left) + (padding * 2),
            height: (bottom - top) + (padding * 2)
        };

    });

    if (clip) {
        await page.screenshot({ path: screenshotPath, clip });
        console.log(`Saved cropped image: ${screenshotPath}`);
    } else {
        console.warn(`Skipping ${id}, elements not found.`);
    }
}




    await browser.close();
})();
