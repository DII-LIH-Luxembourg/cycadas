// Marker selector input --------------------------------------------------------
// A searchable grid of markers; each marker is off, positive or negative.
// Value sent to Shiny: {pos: [...], neg: [...]} in marker order.
// Server messages (update_marker_selector()):
//   markers  array of marker names
//   locked   {marker: "pos" | "neg"} markers fixed by the node lineage
//   flagged  array of markers to mark with a warning dot (not bimodal)
//   clear    true to reset all selections

(function () {
  const binding = new Shiny.InputBinding();

  const rowHtml = (m) => `
    <div class="marker-row" data-marker="${m}" data-state="off">
      <span class="marker-name" title="${m}">${m}</span>
      <span class="marker-toggle">
        <button type="button" class="marker-btn" data-set="neg" title="${m} negative">&minus;</button>
        <button type="button" class="marker-btn" data-set="pos" title="${m} positive">+</button>
      </span>
    </div>`;

  $.extend(binding, {
    find: (scope) => $(scope).find(".marker-selector"),

    initialize: function (el) {
      $(el).data("markers", []);
    },

    getValue: function (el) {
      const pick = (state) =>
        $(el).find(`.marker-row[data-state="${state}"]:not(.locked)`)
          .map((i, r) => r.dataset.marker).get();
      return { pos: pick("pos"), neg: pick("neg") };
    },

    subscribe: function (el, callback) {
      $(el).on("click.markerSelector", ".marker-btn", function () {
        const row = this.closest(".marker-row");
        if (row.classList.contains("locked")) return;
        const set = this.dataset.set;
        row.dataset.state = row.dataset.state === set ? "off" : set;
        binding.refreshCount(el);
        callback();
      });
      $(el).on("click.markerSelector", ".marker-clear", function () {
        $(el).find(".marker-row:not(.locked)").attr("data-state", "off");
        binding.refreshCount(el);
        callback();
      });
      $(el).on("input.markerSelector", ".marker-search", function () {
        const q = this.value.trim().toLowerCase();
        $(el).find(".marker-row").each(function () {
          this.hidden = q !== "" && !this.dataset.marker.toLowerCase().includes(q);
        });
      });
      $(el).on("change.markerSelector", () => callback());
    },

    unsubscribe: (el) => $(el).off(".markerSelector"),

    refreshCount: function (el) {
      const v = binding.getValue(el);
      const n = v.pos.length + v.neg.length;
      $(el).find(".marker-count").text(n ? `${n} selected` : "");
      $(el).find(".marker-clear").prop("hidden", n === 0);
    },

    receiveMessage: function (el, msg) {
      const grid = $(el).find(".marker-grid");
      if (msg.markers) {
        const keep = msg.clear ? {} : Object.fromEntries(
          $(el).find(".marker-row").map((i, r) => [[r.dataset.marker, r.dataset.state]]).get());
        grid.html(msg.markers.map(rowHtml).join(""));
        grid.find(".marker-row").each(function () {
          if (keep[this.dataset.marker]) this.dataset.state = keep[this.dataset.marker];
        });
        $(el).find(".marker-empty").prop("hidden", msg.markers.length > 0);
      } else if (msg.clear) {
        grid.find(".marker-row").attr("data-state", "off");
      }
      if (msg.locked) {
        grid.find(".marker-row").each(function () {
          const lock = msg.locked[this.dataset.marker];
          this.classList.toggle("locked", !!lock);
          if (lock) this.dataset.state = lock;
          else if (this.dataset.state !== "off" && msg.clear) this.dataset.state = "off";
        });
      }
      if (msg.flagged) {
        const flagged = new Set(msg.flagged);
        grid.find(".marker-row").each(function () {
          this.classList.toggle("flagged", flagged.has(this.dataset.marker));
        });
      }
      $(el).find(".marker-search").trigger("input");
      binding.refreshCount(el);
      $(el).trigger("change");
    }
  });

  Shiny.inputBindings.register(binding, "cycadas.markerSelector");
})();

// Tree view controls -----------------------------------------------------------
// Buttons from tree_toolbar() call cycadasTree.run(widgetId, action) on the
// visNetwork widget; readable() is also used for the initial view.
window.cycadasTree = (function () {
  const MIN_SCALE = 0.65, MAX_SCALE = 1, ROOT_ID = 1;
  const anim = { duration: 250, easingFunction: "easeInOutQuad" };

  // Fit the tree, but zoom out no further than MIN_SCALE (large trees, anchored
  // at the root on the left) and zoom in no further than MAX_SCALE (small trees)
  function readable(network, animation) {
    network.fit({ animation: false });
    const scale = network.getScale();
    if (scale < MIN_SCALE) {
      const root = network.getPositions([ROOT_ID])[ROOT_ID];
      const width = network.body.container.clientWidth;
      const position = root ? { x: root.x + width / 2 / MIN_SCALE - 90, y: root.y }
                            : network.getViewPosition();
      network.moveTo({ scale: MIN_SCALE, position: position, animation: animation || false });
    } else if (scale > MAX_SCALE) {
      network.moveTo({ scale: MAX_SCALE, animation: animation || false });
    }
  }

  function pan(network, dx, dy) {
    const step = 150 / network.getScale();
    const p = network.getViewPosition();
    network.moveTo({ position: { x: p.x + dx * step, y: p.y + dy * step }, animation: anim });
  }

  function run(id, action) {
    // visNetwork keeps the vis.js network on an inner element "graph<id>"
    const el = document.getElementById("graph" + id);
    const network = el && el.chart;
    if (!network) return;
    const scale = network.getScale();
    switch (action) {
      case "zoomIn":  network.moveTo({ scale: scale * 1.25, animation: anim }); break;
      case "zoomOut": network.moveTo({ scale: scale / 1.25, animation: anim }); break;
      case "fit":     network.fit({ animation: anim }); break;
      case "reset":   readable(network); break;
      case "focus": {
        const sel = network.getSelectedNodes();
        if (sel.length) network.focus(sel[0], { scale: Math.max(scale, MAX_SCALE), animation: anim });
        break;
      }
      case "left":  pan(network, -1, 0); break;
      case "right": pan(network, 1, 0); break;
      case "up":    pan(network, 0, -1); break;
      case "down":  pan(network, 0, 1); break;
    }
  }

  // SVG of the whole tree with the layout shown on screen ---------------------
  // Edges follow the horizontal cubic Bezier curves drawn by vis.js
  // (roundness 0.5); boxes use the node bounding boxes.
  const esc = (t) => String(t).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

  function buildSvg(network) {
    const nodes = network.body.data.nodes.get();
    const edges = network.body.data.edges.get();
    const pos = network.getPositions();
    const box = {};
    nodes.forEach((n) => { box[n.id] = network.getBoundingBox(n.id); });

    const pad = 20;
    const minX = Math.min(...nodes.map((n) => box[n.id].left)) - pad;
    const minY = Math.min(...nodes.map((n) => box[n.id].top)) - pad;
    const width = Math.max(...nodes.map((n) => box[n.id].right)) - minX + pad;
    const height = Math.max(...nodes.map((n) => box[n.id].bottom)) - minY + pad;
    const X = (x) => (x - minX).toFixed(1), Y = (y) => (y - minY).toFixed(1);

    const paths = edges.filter((e) => pos[e.from] && pos[e.to]).map((e) => {
      const a = pos[e.from], b = pos[e.to], dx = (b.x - a.x) * 0.5;
      return `<path d="M${X(a.x)},${Y(a.y)} C${X(a.x + dx)},${Y(a.y)} ${X(b.x - dx)},${Y(b.y)} ${X(b.x)},${Y(b.y)}"/>`;
    });

    const boxes = nodes.map((n) => {
      const bb = box[n.id];
      const fill = n["color.background"] || (n.color && n.color.background) || "#2c7be5";
      const text = n["font.color"] || "#ffffff";
      const cy = (bb.top + bb.bottom) / 2;
      return `<g><rect x="${X(bb.left)}" y="${Y(bb.top)}" width="${(bb.right - bb.left).toFixed(1)}" ` +
        `height="${(bb.bottom - bb.top).toFixed(1)}" rx="6" fill="${fill}"/>` +
        `<text x="${X(pos[n.id].x)}" y="${Y(cy)}" fill="${text}">${esc(n.label)}</text></g>`;
    });

    return `<svg xmlns="http://www.w3.org/2000/svg" width="${width.toFixed(0)}" height="${height.toFixed(0)}" ` +
      `viewBox="0 0 ${width.toFixed(1)} ${height.toFixed(1)}">` +
      `<rect width="100%" height="100%" fill="#ffffff"/>` +
      `<g fill="none" stroke="#ced4da" stroke-width="1.5">${paths.join("")}</g>` +
      `<g font-family="system-ui, -apple-system, 'Segoe UI', Roboto, Arial, sans-serif" font-size="15" ` +
      `text-anchor="middle" dominant-baseline="central">${boxes.join("")}</g></svg>`;
  }

  function download(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  function exportTree(id, format) {
    const el = document.getElementById("graph" + id);
    const network = el && el.chart;
    if (!network) return;
    const svg = buildSvg(network);
    const name = "cycadas_tree_" + new Date().toISOString().slice(0, 10);
    if (format === "svg") {
      download(new Blob([svg], { type: "image/svg+xml" }), name + ".svg");
      return;
    }
    // PNG at 3x resolution, rendered from the SVG
    const img = new Image();
    img.onload = function () {
      const canvas = document.createElement("canvas");
      canvas.width = img.width * 3; canvas.height = img.height * 3;
      const ctx = canvas.getContext("2d");
      ctx.scale(3, 3);
      ctx.drawImage(img, 0, 0);
      canvas.toBlob((blob) => download(blob, name + ".png"), "image/png");
    };
    img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(svg);
  }

  return { readable: readable, run: run, buildSvg: buildSvg, exportTree: exportTree };
})();
