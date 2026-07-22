/* ==========================================================================
   aucmat pkgdown extras — GLightbox, reading bar, back-to-top
   Uses GLightbox (same as Quarto) for click-to-zoom with smooth animation.
   ========================================================================== */

(function() {

  /* ========================================================================
     Create DOM elements for progress bar and back-to-top
     ======================================================================== */

  var rpb = document.createElement('div');
  rpb.id = 'rpb';
  document.body.appendChild(rpb);

  var btt = document.createElement('button');
  btt.id = 'btt';
  btt.title = 'Back to top';
  btt.setAttribute('aria-label', 'Back to top');
  btt.innerHTML = '↑';
  document.body.appendChild(btt);

  /* ── Progress bar ──────────────────────────────── */
  window.addEventListener('scroll', function() {
    var s = document.documentElement;
    rpb.style.width = Math.min(
      (s.scrollTop / (s.scrollHeight - s.clientHeight)) * 100, 100
    ) + '%';
  }, { passive: true });

  /* ── Back to top ───────────────────────────────── */
  window.addEventListener('scroll', function() {
    btt.classList.toggle('show', window.scrollY > 480);
  }, { passive: true });
  btt.addEventListener('click', function() {
    window.scrollTo({ top: 0, behavior: 'smooth' });
  });

  /* ========================================================================
     GLightbox — load from CDN, init on all content images
     ======================================================================== */

  function loadGLightbox(cb) {
    /* CSS */
    var css = document.createElement('link');
    css.rel = 'stylesheet';
    css.href = 'https://cdn.jsdelivr.net/npm/glightbox/dist/css/glightbox.min.css';
    document.head.appendChild(css);

    /* JS */
    var js = document.createElement('script');
    js.src = 'https://cdn.jsdelivr.net/npm/glightbox/dist/js/glightbox.min.js';
    js.onload = cb;
    document.head.appendChild(js);
  }

  function initLightbox() {
    /* Wrap qualifying images in <a class="lightbox"> links so GLightbox
       can bind to them.  Skip navbar icons, tiny images, and SVGs. */
    var images = document.querySelectorAll(
      'main img[src], article img[src], .figure img[src], ' +
      '.section img[src], .content img[src]'
    );

    images.forEach(function(img) {
      if (img.dataset.lbx) return;
      if (img.closest && (img.closest('nav') || img.closest('.navbar') ||
          img.closest('.glightbox'))) return;
      if (!img.src || img.src.indexOf('data:image/svg') === 0) return;
      if (img.naturalWidth < 60 && img.naturalHeight < 60) return;

      img.dataset.lbx = '1';

      /* Create wrapper <a> if not already inside one */
      var parent = img.parentNode;
      if (parent.tagName === 'A') {
        parent.classList.add('lightbox');
        return;
      }

      var a = document.createElement('a');
      a.href = img.src;
      a.className = 'lightbox';
      a.setAttribute('data-gallery', 'aucmat-gallery');
      if (img.alt) a.setAttribute('data-title', img.alt);
      parent.insertBefore(a, img);
      a.appendChild(img);
    });

    /* Init GLightbox with Quarto-like settings */
    if (typeof GLightbox !== 'undefined') {
      window._aucmatLightbox = GLightbox({
        closeEffect: 'zoom',
        descPosition: 'bottom',
        loop: false,
        openEffect: 'zoom',
        selector: '.lightbox'
      });
    }
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', function() {
      loadGLightbox(initLightbox);
    });
  } else {
    loadGLightbox(initLightbox);
  }

})();
