/* ==========================================================================
   aucmat pkgdown extras — lightbox, reading bar, back-to-top
   Runs after DOM is ready.
   ========================================================================== */

(function() {

  /* ========================================================================
     Create DOM elements for lightbox, progress bar, back-to-top
     ======================================================================== */

  /* ── Reading progress bar ──────────────────────── */
  var rpb = document.createElement('div');
  rpb.id = 'rpb';
  document.body.appendChild(rpb);

  /* ── Back to top button ────────────────────────── */
  var btt = document.createElement('button');
  btt.id = 'btt';
  btt.title = 'Back to top';
  btt.setAttribute('aria-label', 'Back to top');
  btt.innerHTML = '↑';  /* ↑ */
  document.body.appendChild(btt);

  /* ── Lightbox overlay ──────────────────────────── */
  var overlay = document.createElement('div');
  overlay.id = 'lbx-overlay';
  var lbxImg = document.createElement('img');
  lbxImg.id = 'lbx-img';
  lbxImg.setAttribute('draggable', 'false');
  overlay.appendChild(lbxImg);

  var hint = document.createElement('div');
  hint.id = 'lbx-hint';
  hint.textContent = 'Click to zoom in • Scroll to adjust • Esc to close';

  /* ── Control buttons ───────────────────────────── */
  var btnClose = document.createElement('button');
  btnClose.id = 'lbx-close';
  btnClose.title = 'Close (Esc)';
  btnClose.innerHTML = '×';  /* × */

  var btnIn = document.createElement('button');
  btnIn.id = 'lbx-zoom-in';
  btnIn.title = 'Zoom in (+)';
  btnIn.innerHTML = '+';

  var btnOut = document.createElement('button');
  btnOut.id = 'lbx-zoom-out';
  btnOut.title = 'Zoom out (−)';
  btnOut.innerHTML = '−';  /* − */

  var btnFit = document.createElement('button');
  btnFit.id = 'lbx-fit';
  btnFit.title = 'Fit to screen (0)';
  btnFit.innerHTML = '⦙';  /* ⤡ */

  var btnFull = document.createElement('button');
  btnFull.id = 'lbx-full';
  btnFull.title = 'Actual size (1)';
  btnFull.innerHTML = '1:1';

  document.body.appendChild(overlay);
  document.body.appendChild(hint);
  document.body.appendChild(btnClose);
  document.body.appendChild(btnIn);
  document.body.appendChild(btnOut);
  document.body.appendChild(btnFit);
  document.body.appendChild(btnFull);

  /* ========================================================================
     Logic
     ======================================================================== */

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

  /* ── Lightbox state ────────────────────────────── */
  var scale = 1, panX = 0, panY = 0,
      fitScale = 1, naturalW = 0, naturalH = 0,
      isPanning = false, lastX = 0, lastY = 0,
      MIN_SCALE = 0.2, MAX_SCALE = 10, ZOOM_STEP = 1.4;

  function showHint(msg) {
    hint.textContent = msg || 'Click to zoom in • Scroll to adjust • Esc to close';
    hint.classList.add('show');
    clearTimeout(hint._t);
    hint._t = setTimeout(function() { hint.classList.remove('show'); }, 2500);
  }

  function applyTransform(s, x, y) {
    var tx = (x || 0), ty = (y || 0);
    lbxImg.style.transform = 'translate(' + tx + 'px,' + ty + 'px) scale(' + s + ')';
    var zoomed = s > fitScale * 1.02;
    lbxImg.classList.toggle('zoomed', zoomed);
  }

  function fitToScreen() {
    var vw = window.innerWidth * 0.9, vh = window.innerHeight * 0.85;
    fitScale = Math.min(vw / naturalW, vh / naturalH, 2);
    scale = fitScale; panX = 0; panY = 0;
    applyTransform(fitScale, 0, 0);
    lbxImg.style.cursor = 'zoom-in';
    showHint('Click to see full size');
  }

  function zoomToFull() {
    scale = 1; panX = 0; panY = 0;
    applyTransform(1, 0, 0);
    lbxImg.style.cursor = 'grab';
    showHint('Scroll to zoom • Drag to pan • Click to fit screen');
  }

  function open(src) {
    lbxImg.src = src;
    var tmp = new Image();
    tmp.onload = function() {
      naturalW = tmp.naturalWidth; naturalH = tmp.naturalHeight;
      fitToScreen();
    };
    tmp.src = src;
    if (tmp.complete) {
      naturalW = tmp.naturalWidth; naturalH = tmp.naturalHeight;
      fitToScreen();
    }
    overlay.classList.add('active');
    document.body.style.overflow = 'hidden';
    showHint();
  }

  function close() {
    overlay.classList.remove('active');
    document.body.style.overflow = '';
    scale = 1; panX = 0; panY = 0;
    lbxImg.style.transform = '';
    lbxImg.style.cursor = 'zoom-in';
    lbxImg.classList.remove('zoomed', 'panning');
    hint.classList.remove('show');
  }

  function clampPan(s, x, y) {
    var imgW = naturalW * s, imgH = naturalH * s,
        vw = window.innerWidth, vh = window.innerHeight,
        maxX = Math.max(0, (imgW - vw) / 2),
        maxY = Math.max(0, (imgH - vh) / 2);
    return [
      Math.max(-maxX, Math.min(maxX, x || 0)),
      Math.max(-maxY, Math.min(maxY, y || 0))
    ];
  }

  /* ── Click image: toggle fit ↔ full ───────────── */
  lbxImg.addEventListener('click', function(e) {
    if (isPanning) return;
    e.stopPropagation();
    if (scale > fitScale * 1.02) { fitToScreen(); }
    else { zoomToFull(); }
  });

  /* ── Double-click: straight to actual size ─────── */
  lbxImg.addEventListener('dblclick', function(e) {
    e.stopPropagation(); e.preventDefault();
    zoomToFull();
  });

  /* ── Click background: close ───────────────────── */
  overlay.addEventListener('click', function(e) {
    if (e.target === overlay) close();
  });

  /* ── Scroll wheel: zoom centered on cursor ─────── */
  overlay.addEventListener('wheel', function(e) {
    if (!overlay.classList.contains('active')) return;
    e.preventDefault();
    var rect = lbxImg.getBoundingClientRect(),
        mx = e.clientX - rect.left - rect.width / 2,
        my = e.clientY - rect.top - rect.height / 2,
        oldScale = scale;

    if (e.deltaY < 0) scale = Math.min(MAX_SCALE, scale * ZOOM_STEP);
    else scale = Math.max(MIN_SCALE, scale / ZOOM_STEP);

    var ratio = scale / oldScale;
    var clamped = clampPan(scale, panX * ratio + mx * (1 - ratio),
                                  panY * ratio + my * (1 - ratio));
    panX = clamped[0]; panY = clamped[1];
    applyTransform(scale, panX, panY);
    lbxImg.style.cursor = scale > fitScale * 1.02 ? 'grab' : 'zoom-in';
  }, { passive: false });

  /* ── Drag to pan ───────────────────────────────── */
  lbxImg.addEventListener('mousedown', function(e) {
    if (scale <= fitScale * 1.02) return;
    e.preventDefault();
    isPanning = true; lastX = e.clientX; lastY = e.clientY;
    lbxImg.classList.add('panning');
  });
  window.addEventListener('mousemove', function(e) {
    if (!isPanning) return;
    var dx = e.clientX - lastX, dy = e.clientY - lastY;
    lastX = e.clientX; lastY = e.clientY;
    var clamped = clampPan(scale, panX + dx, panY + dy);
    panX = clamped[0]; panY = clamped[1];
    applyTransform(scale, panX, panY);
  });
  window.addEventListener('mouseup', function() {
    if (isPanning) { isPanning = false; lbxImg.classList.remove('panning'); }
  });

  /* ── Touch pinch zoom ──────────────────────────── */
  var pinchDist0 = 0, pinchScale0 = 1;
  overlay.addEventListener('touchstart', function(e) {
    if (e.touches.length === 2) {
      pinchDist0 = Math.hypot(
        e.touches[0].clientX - e.touches[1].clientX,
        e.touches[0].clientY - e.touches[1].clientY
      );
      pinchScale0 = scale;
    }
  }, { passive: false });
  overlay.addEventListener('touchmove', function(e) {
    if (e.touches.length === 2 && overlay.classList.contains('active')) {
      e.preventDefault();
      var dist = Math.hypot(
        e.touches[0].clientX - e.touches[1].clientX,
        e.touches[0].clientY - e.touches[1].clientY
      );
      var s = Math.min(MAX_SCALE, Math.max(MIN_SCALE, pinchScale0 * (dist / pinchDist0)));
      applyTransform(s, panX, panY);
      scale = s;
    }
  }, { passive: false });

  /* ── Button handlers ───────────────────────────── */
  btnClose.addEventListener('click', close);

  btnIn.addEventListener('click', function() {
    scale = Math.min(MAX_SCALE, scale * ZOOM_STEP);
    var c = clampPan(scale, panX, panY);
    panX = c[0]; panY = c[1];
    applyTransform(scale, panX, panY);
    lbxImg.style.cursor = 'grab';
  });

  btnOut.addEventListener('click', function() {
    scale = Math.max(MIN_SCALE, scale / ZOOM_STEP);
    var c = clampPan(scale, panX, panY);
    panX = c[0]; panY = c[1];
    applyTransform(scale, panX, panY);
    if (scale <= fitScale * 1.02) lbxImg.style.cursor = 'zoom-in';
  });

  btnFit.addEventListener('click', function() {
    fitToScreen(); lbxImg.style.cursor = 'zoom-in';
  });

  btnFull.addEventListener('click', function() {
    zoomToFull();
  });

  /* ── Keyboard shortcuts ────────────────────────── */
  document.addEventListener('keydown', function(e) {
    if (!overlay.classList.contains('active')) return;
    if (e.key === 'Escape') { close(); return; }
    if (e.key === '+' || e.key === '=') {
      scale = Math.min(MAX_SCALE, scale * ZOOM_STEP);
      applyTransform(scale, panX, panY);
      lbxImg.style.cursor = 'grab';
    }
    if (e.key === '-') {
      scale = Math.max(MIN_SCALE, scale / ZOOM_STEP);
      applyTransform(scale, panX, panY);
    }
    if (e.key === '0') { fitToScreen(); lbxImg.style.cursor = 'zoom-in'; }
    if (e.key === '1') { zoomToFull(); }
  });

  /* ── Hook: click any content image opens lightbox ─ */
  function hookImages() {
    var images = document.querySelectorAll('main img[src], article img[src]');
    images.forEach(function(img) {
      if (img.id === 'lbx-img') return;
      if (img.closest && (img.closest('#lbx-overlay') ||
          img.closest('nav') || img.closest('.navbar'))) return;
      if (!img.src || img.src.indexOf('data:image/svg') === 0) return;
      if (img.naturalWidth < 60 && img.naturalHeight < 60) return;
      if (img.dataset.lbxHooked) return;
      img.dataset.lbxHooked = '1';
      img.style.cursor = 'zoom-in';
      img.addEventListener('click', function(e) {
        e.preventDefault(); open(img.src);
      });
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', hookImages);
  } else {
    hookImages();
  }

  /* Re-hook after pkgdown page transitions (bootstrap nav) */
  if (typeof window.$ !== 'undefined') {
    window.$(document).on('shown.bs.tab', hookImages);
  }

})();
