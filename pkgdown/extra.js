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
    btt.innerHTML = '↑';
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

    /* ── Wrap images in .lightbox <a> tags ───────── */
    function wrapImages() {
      if (wrapImages._done) return;
      var imgs = document.querySelectorAll('img[src]');
      for (var i = 0; i < imgs.length; i++) {
        var img = imgs[i];
        if (img.dataset.lbx) continue;
        var p = img.closest && img.closest('nav,.navbar,.glightbox');
        if (p) continue;
        if (!img.src || /^data:/.test(img.src)) continue;
        if (img.naturalWidth < 50 && img.naturalHeight < 50) continue;

        img.dataset.lbx = '1';
        var parent = img.parentNode;
        if (!parent) continue;

        if (parent.tagName === 'A') {
          parent.classList.add('lightbox');
          if (!parent.hasAttribute('data-gallery')) parent.setAttribute('data-gallery', 'aucmat');
          continue;
        }

        var a = document.createElement('a');
        a.href = img.src;
        a.className = 'lightbox';
        a.setAttribute('data-gallery', 'aucmat');
        if (img.alt) a.setAttribute('data-title', img.alt);
        parent.insertBefore(a, img);
        a.appendChild(img);
      }
      wrapImages._done = true;
    }

    /* ── Init GLightbox ──────────────────────────── */
    var _initialized = false;
    function initGL() {
      if (_initialized) return;
      if (typeof GLightbox === 'undefined') return;
      _initialized = true;
      wrapImages();
      window._lbx = GLightbox({
        closeEffect: 'zoom',
        descPosition: 'bottom',
        loop: false,
        openEffect: 'zoom',
        selector: '.lightbox'
      });
    }

    /* Load GLightbox JS from CDN */
    var js = document.createElement('script');
    js.src = 'https://cdn.jsdelivr.net/npm/glightbox/dist/js/glightbox.min.js';
    js.onload = initGL;
    document.head.appendChild(js);

    /* Fallback: poll if onload doesn't fire */
    var pollTries = 0;
    function poll() {
      if (_initialized) return;
      if (typeof GLightbox !== 'undefined') { initGL(); return; }
      if (++pollTries < 60) setTimeout(poll, 100);
    }
    setTimeout(poll, 1000);

  }); /* end whenReady */

})();
