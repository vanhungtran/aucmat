// KaTeX auto-render for Pandoc-style math spans
// Handles <span class="math inline"> and <span class="math display">
// Reference: https://github.com/jgm/pandoc/blob/29fa97ab96b8e2d62d48326e1b949a71dc41f47a/src/Text/Pandoc/Writers/HTML.hs#L332-L345
(function () {
  function renderMath() {
    if (typeof katex === 'undefined') {
      // KaTeX CDN not loaded yet; retry once
      setTimeout(renderMath, 100);
      return;
    }
    var mathElements = document.getElementsByClassName("math");
    for (var i = 0; i < mathElements.length; i++) {
      var el = mathElements[i];
      // Skip already-rendered elements
      if (el.querySelector('.katex') || el.querySelector('.katex-error')) continue;
      if (el.tagName !== "SPAN") continue;
      var texText = el.firstChild;
      if (!texText || !texText.data) continue;
      try {
        katex.render(texText.data, el, {
          displayMode: el.classList.contains("display"),
          throwOnError: false,
          fleqn: false
        });
      } catch (e) {
        // Leave raw LaTeX visible on render failure
        console.warn('KaTeX render failed:', e.message);
      }
    }
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', renderMath);
  } else {
    renderMath();
  }
})();
