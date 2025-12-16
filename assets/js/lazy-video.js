(function () {
  function loadVideo(v) {
    if (v.dataset.loaded) return;

    const src = v.dataset.src;
    if (!src) return;

    v.src = src;
    v.dataset.loaded = "1";
    v.load();
  }

  function inOrNearViewport(el, preloadPx) {
    const rect = el.getBoundingClientRect();
    const vpH = window.innerHeight || document.documentElement.clientHeight;
    const vpW = window.innerWidth || document.documentElement.clientWidth;

    return (
      rect.bottom >= -preloadPx &&
      rect.right >= 0 &&
      rect.top <= vpH + preloadPx &&
      rect.left <= vpW
    );
  }

  function init() {
    const videos = document.querySelectorAll("video.lazy-video");
    if (!videos.length) return;

    const preloadPx = 200;

    // Always load anything already in/near view (fixes cached fast-load cases)
    videos.forEach((v) => {
      if (inOrNearViewport(v, preloadPx)) loadVideo(v);
    });

    if (!("IntersectionObserver" in window)) {
      videos.forEach(loadVideo);
      return;
    }

    const io = new IntersectionObserver(
      (entries) => {
        for (const entry of entries) {
          if (entry.isIntersecting) {
            loadVideo(entry.target);
            io.unobserve(entry.target);
          }
        }
      },
      { rootMargin: preloadPx + "px 0px", threshold: 0.01 }
    );

    videos.forEach((v) => {
      if (!v.dataset.loaded) io.observe(v);
    });

    // Extra safety: sometimes layout shifts after DOMContentLoaded
    // (fonts/images/MathJax). Recheck once after full load.
    window.addEventListener(
      "load",
      () => {
        videos.forEach((v) => {
          if (!v.dataset.loaded && inOrNearViewport(v, preloadPx)) {
            loadVideo(v);
            io.unobserve(v);
          }
        });
      },
      { once: true }
    );
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
