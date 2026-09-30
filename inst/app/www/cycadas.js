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
