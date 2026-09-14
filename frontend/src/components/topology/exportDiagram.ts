/** Export SVG diagram as high-resolution PNG or JPEG for publication figures. */

const EXPORT_STYLES = `
  .membrane-label { fill: #64748b; font-size: 10px; font-weight: 700; letter-spacing: .12em; font-family: Arial, Helvetica, sans-serif; }
  .helix-label-text { fill: #ffffff; font-size: 13px; font-weight: 800; font-family: Arial, Helvetica, sans-serif; }
  .helix-sub-text { fill: #f8fafc; font-size: 8px; font-weight: 700; font-family: Arial, Helvetica, sans-serif; }
  .loop-text-label { fill: #1e293b; font-size: 11px; font-weight: 700; font-family: Arial, Helvetica, sans-serif; }
  .short-helix-text { fill: #ffffff; font-size: 9px; font-weight: 800; font-family: Arial, Helvetica, sans-serif; }
  .terminal-text { fill: #0f766e; font-size: 19px; font-weight: 800; font-family: Georgia, serif; }
  .domain-track-label { fill: #0f766e; font-size: 12px; font-weight: 800; font-family: Arial, Helvetica, sans-serif; }
  .domain-bar-text { fill: #ffffff; font-size: 9px; font-weight: 800; font-family: Arial, Helvetica, sans-serif; }
  .panel-letter { fill: #1e293b; font-size: 22px; font-weight: 800; font-family: Arial, Helvetica, sans-serif; }
`;

export async function exportSvgAsImage(
  svgElement: SVGSVGElement,
  format: 'png' | 'jpeg',
  filename: string,
  scale = 3,
): Promise<void> {
  const clone = svgElement.cloneNode(true) as SVGSVGElement;

  const styleEl = document.createElementNS('http://www.w3.org/2000/svg', 'style');
  styleEl.textContent = EXPORT_STYLES;
  clone.insertBefore(styleEl, clone.firstChild);

  const vb = svgElement.viewBox.baseVal;
  const vbWidth = vb.width || svgElement.clientWidth || 940;
  const vbHeight = vb.height || svgElement.clientHeight || 440;

  clone.setAttribute('width', String(vbWidth));
  clone.setAttribute('height', String(vbHeight));

  const svgData = new XMLSerializer().serializeToString(clone);
  const svgBlob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(svgBlob);

  try {
    await new Promise<void>((resolve, reject) => {
      const img = new Image();
      img.onload = () => {
        const canvas = document.createElement('canvas');
        canvas.width = vbWidth * scale;
        canvas.height = vbHeight * scale;
        const ctx = canvas.getContext('2d');
        if (!ctx) {
          reject(new Error('Canvas not supported'));
          return;
        }

        if (format === 'jpeg') {
          ctx.fillStyle = '#ffffff';
          ctx.fillRect(0, 0, canvas.width, canvas.height);
        }

        ctx.drawImage(img, 0, 0, canvas.width, canvas.height);

        canvas.toBlob(
          (blob) => {
            if (!blob) {
              reject(new Error('Export failed'));
              return;
            }
            const link = document.createElement('a');
            link.href = URL.createObjectURL(blob);
            link.download = filename;
            link.click();
            URL.revokeObjectURL(link.href);
            resolve();
          },
          format === 'jpeg' ? 'image/jpeg' : 'image/png',
          0.95,
        );
      };
      img.onerror = () => reject(new Error('Failed to render SVG'));
      img.src = url;
    });
  } finally {
    URL.revokeObjectURL(url);
  }
}
