// inline.js - Inlines dist/cheetah.js directly into demo/index.html and demo/3d.html to create self-contained pages

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const bundledJsPath = path.join(__dirname, '../dist/cheetah.js');

function inlineFile(inputHtmlName, outputHtmlName) {
    const demoHtmlPath = path.join(__dirname, `../demo/${inputHtmlName}`);
    const outputHtmlPath = path.join(__dirname, `../docs/${outputHtmlName}`);
    
    console.log(`Generating self-contained page ${outputHtmlName}...`);
    
    // Ensure docs directory exists
    const docsDir = path.dirname(outputHtmlPath);
    if (!fs.existsSync(docsDir)) {
        fs.mkdirSync(docsDir, { recursive: true });
    }
    
    if (!fs.existsSync(demoHtmlPath)) {
        console.error(`Demo HTML not found at: ${demoHtmlPath}`);
        process.exit(1);
    }
    if (!fs.existsSync(bundledJsPath)) {
        console.error(`Bundled JS not found at: ${bundledJsPath}. Run 'npm run build' first.`);
        process.exit(1);
    }

    let html = fs.readFileSync(demoHtmlPath, 'utf8');
    const js = fs.readFileSync(bundledJsPath, 'utf8');

    // Replace the external script tag with inlined script content
    const scriptRegex = /<script\s+src="\.\.\/dist\/cheetah\.js"><\/script>/i;
    
    if (!scriptRegex.test(html)) {
        console.error(`Could not find reference script tag in ${inputHtmlName}`);
        process.exit(1);
    }

    const inlinedScript = `<script>\n${js}\n</script>`;
    html = html.replace(scriptRegex, inlinedScript);

    fs.writeFileSync(outputHtmlPath, html, 'utf8');
    console.log(`Successfully generated self-contained HTML page at: ${outputHtmlPath}`);
}

inlineFile('index.html', 'index.html');
inlineFile('3d.html', '3d.html');
