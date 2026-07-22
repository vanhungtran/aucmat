/* ==========================================================================
   aucmat pkgdown extras — GLightbox zoom, reading bar, back-to-top
   ========================================================================== */
(function() {

  /* ── All DOM operations deferred until body exists ─ */
  function whenReady(fn) {
    if (document.body) { fn(); return; }
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', fn);
    } else {
      setTimeout(function() { whenReady(fn); }, 10);
    }
  }

  whenReady(function() {

    /* ── Progress bar + back-to-top ──────────────── */
    var rpb = document.createElement('div'), btt = document.createElement('button');
    rpb.id = 'rpb'; btt.id = 'btt';
    btt.title = 'Back to top'; btt.setAttribute('aria-label', 'Back to top');
    btt.innerHTML = '\u2191';
    document.body.appendChild(rpb);
    document.body.appendChild(btt);

    window.addEventListener('scroll', function() {
      var s = document.documentElement;
      rpb.style.width = Math.min((s.scrollTop / (s.scrollHeight - s.clientHeight)) * 100, 100) + '%';
      btt.classList.toggle('show', window.scrollY > 480);
    }, { passive: true });
    btt.addEventListener('click', function() { window.scrollTo({ top: 0, behavior: 'smooth' }); });

    /* ── Load GLightbox CSS ──────────────────────── */
    var css = document.createElement('link');
    css.rel = 'stylesheet';
    css.href = 'https://cdn.jsdelivr.net/npm/glightbox/dist/css/glightbox.min.css';
    document.head.appendChild(css);

    /* ── Process a single image: wrap in .lightbox <a> tag ─ */
    function processImage(img) {
      if (img.dataset.lbx) return;
      if (!img.src || /^data:/.test(img.src)) return;
      // skip tiny icons / inline SVGs
      if (img.naturalWidth > 0 && img.naturalHeight > 0 &&
          img.naturalWidth < 50 && img.naturalHeight < 50) return;

      // skip nav, navbar, and already-lightbox containers
      if (img.closest && img.closest('nav,.navbar,.glightbox,.goverlay')) return;

      img.dataset.lbx = '1';
      var parent = img.parentNode;
      if (!parent) return;

      if (parent.tagName === 'A') {
        parent.classList.add('lightbox');
        if (!parent.hasAttribute('data-gallery')) parent.setAttribute('data-gallery', 'aucmat');
        return;
      }

      var a = document.createElement('a');
      a.href = img.src;
      a.className = 'lightbox';
      a.setAttribute('data-gallery', 'aucmat');
      if (img.alt) a.setAttribute('data-title', img.alt);
      parent.insertBefore(a, img);
      a.appendChild(img);
    }

    /* ── Scan all existing images ──────────────── */
    function scanAllImages() {
      var imgs = document.querySelectorAll('img[src]');
      for (var i = 0; i < imgs.length; i++) {
        processImage(imgs[i]);
      }
    }

    /* ── MutationObserver: catch images added later ─ */
    function watchForImages() {
      if (watchForImages._active) return;
      watchForImages._active = true;
      var observer = new MutationObserver(function(mutations) {
        for (var m = 0; m < mutations.length; m++) {
          var nodes = mutations[m].addedNodes;
          for (var n = 0; n < nodes.length; n++) {
            var node = nodes[n];
            if (node.nodeType !== 1) continue;
            if (node.tagName === 'IMG') { processImage(node); }
            if (node.querySelectorAll) {
              var imgs = node.querySelectorAll('img[src]');
              for (var k = 0; k < imgs.length; k++) processImage(imgs[k]);
            }
          }
        }
      });
      observer.observe(document.body || document.documentElement, {
        childList: true, subtree: true
      });
    }

    /* ── Init GLightbox ──────────────────────────── */
    var _initialized = false;
    function initGL() {
      if (_initialized) return;
      if (typeof GLightbox === 'undefined') return;
      _initialized = true;

      scanAllImages();
      watchForImages();

      window._lbx = GLightbox({
        closeEffect: 'zoom',
        descPosition: 'bottom',
        loop: false,
        openEffect: 'zoom',
        selector: '.lightbox'
      });
    }

    /* ── Load GLightbox JS from CDN ────────────── */
    var js = document.createElement('script');
    js.src = 'https://cdn.jsdelivr.net/npm/glightbox/dist/js/glightbox.min.js';
    js.onload = initGL;
    document.head.appendChild(js);

    /* Fallback: poll if CDN onload doesn't fire */
    var pollTries = 0;
    function poll() {
      if (_initialized) return;
      if (typeof GLightbox !== 'undefined') { initGL(); return; }
      if (++pollTries < 60) setTimeout(poll, 100);
    }
    setTimeout(poll, 1000);

    /* Safety net: re-scan after full page load (catches late-rendered plots) */
    window.addEventListener('load', function() {
      setTimeout(function() { scanAllImages(); }, 300);
    });

  }); /* end whenReady */

})();
