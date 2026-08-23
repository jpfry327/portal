/* ChordStats — self-contained statistics for the chordoma atlas static site.
 *
 * Implements: Pearson/Spearman correlation (t-approx p), Kruskal-Wallis
 * (tie-corrected, chi2 p), Kaplan-Meier curves, two-group log-rank test,
 * single-covariate Cox PH with Efron tie handling (matches lifelines),
 * and maximally-selected-log-rank optimal cutpoint.
 *
 * No dependencies. Attaches `ChordStats` to globalThis (browser + jsc).
 */
(function () {
  "use strict";

  // ---------------- special functions ----------------
  function logGamma(x) { // Lanczos, g=7
    const c = [676.5203681218851, -1259.1392167224028, 771.32342877765313,
      -176.61502916214059, 12.507343278686905, -0.13857109526572012,
      9.9843695780195716e-6, 1.5056327351493116e-7];
    if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - logGamma(1 - x);
    x -= 1;
    let a = 0.99999999999980993;
    const t = x + 7.5;
    for (let i = 0; i < 8; i++) a += c[i] / (x + i + 1);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }

  function gammPSeries(a, x) { // regularized lower incomplete gamma, series
    let sum = 1 / a, term = sum;
    for (let n = 1; n < 500; n++) {
      term *= x / (a + n);
      sum += term;
      if (Math.abs(term) < Math.abs(sum) * 1e-15) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - logGamma(a));
  }

  function gammQCF(a, x) { // regularized upper incomplete gamma, continued fraction
    let b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
    for (let i = 1; i < 500; i++) {
      const an = -i * (i - a);
      b += 2;
      d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 1e-15) break;
    }
    return h * Math.exp(-x + a * Math.log(x) - logGamma(a));
  }

  function chi2Sf(x, k) { // P(X > x) for chi-square k df
    if (x <= 0) return 1;
    const a = k / 2, xx = x / 2;
    return xx < a + 1 ? 1 - gammPSeries(a, xx) : gammQCF(a, xx);
  }

  function betacf(a, b, x) { // NR continued fraction for incomplete beta
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < 1e-300) d = 1e-300;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 300; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = 1 + aa / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = 1 + aa / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 1e-14) break;
    }
    return h;
  }

  function ibeta(a, b, x) { // regularized incomplete beta I_x(a,b)
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const bt = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) +
      a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a
      : 1 - bt * betacf(b, a, 1 - x) / b;
  }

  function tSf2(t, df) { // two-sided p for t statistic
    return ibeta(df / 2, 0.5, df / (df + t * t));
  }

  function normSf(z) { // one-sided upper tail of standard normal
    return 0.5 * (z < 0 ? 1 + gammPSeries(0.5, z * z / 2) : gammQCF(0.5, z * z / 2));
  }

  // ---------------- basic stats ----------------
  const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
  function median(a) {          // copies before sorting; does not mutate caller's array
    if (!a.length) return NaN;
    const s = [...a].sort((x, y) => x - y);
    const h = s.length >> 1;
    return s.length % 2 ? s[h] : (s[h - 1] + s[h]) / 2;
  }
  function sd(a) {
    const m = mean(a);
    return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1));
  }

  function ranks(a) { // average ranks for ties
    const idx = a.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]);
    const r = new Array(a.length);
    let i = 0;
    while (i < idx.length) {
      let j = i;
      while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
      const avg = (i + j + 2) / 2;
      for (let k = i; k <= j; k++) r[idx[k][1]] = avg;
      i = j + 1;
    }
    return r;
  }

  function pearson(x, y) {
    const n = x.length, mx = mean(x), my = mean(y);
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < n; i++) {
      sxy += (x[i] - mx) * (y[i] - my);
      sxx += (x[i] - mx) ** 2;
      syy += (y[i] - my) ** 2;
    }
    const r = sxy / Math.sqrt(sxx * syy);
    const t = r * Math.sqrt((n - 2) / (1 - r * r));
    return { r, p: tSf2(Math.abs(t), n - 2), n };
  }

  function spearman(x, y) {
    const res = pearson(ranks(x), ranks(y));
    return { r: res.r, p: res.p, n: x.length };
  }

  function kruskal(groups) { // groups: array of numeric arrays
    const all = [].concat(...groups);
    const N = all.length;
    const r = ranks(all);
    let H = 0, off = 0;
    for (const g of groups) {
      const rg = r.slice(off, off + g.length);
      H += g.length * (mean(rg) - (N + 1) / 2) ** 2;
      off += g.length;
    }
    H *= 12 / (N * (N + 1));
    // tie correction
    const sorted = [...all].sort((a, b) => a - b);
    let tsum = 0, i = 0;
    while (i < N) {
      let j = i;
      while (j + 1 < N && sorted[j + 1] === sorted[i]) j++;
      const t = j - i + 1;
      tsum += t ** 3 - t;
      i = j + 1;
    }
    H /= 1 - tsum / (N ** 3 - N);
    return { H, p: chi2Sf(H, groups.length - 1) };
  }

  // ---------------- survival ----------------
  // Kaplan-Meier estimate. Matches R survfit/ggsurvplot rendering:
  //   - step function drops at event times,
  //   - curve extends horizontally to the LAST observation (censored or not),
  //   - censor marks ("+" in ggsurvplot) returned at censoring times with
  //     the current survival value.
  // Returns { points: [{t,s}], censors: [{t,s}] }.
  function kmCurve(time, event) {
    const idx = time.map((t, i) => i).sort((a, b) => time[a] - time[b]);
    let atRisk = time.length, s = 1;
    const points = [{ t: 0, s: 1 }];
    const censors = [];
    let i = 0, lastT = 0;
    while (i < idx.length) {
      const t = time[idx[i]];
      let d = 0, c = 0;
      while (i < idx.length && time[idx[i]] === t) {
        if (event[idx[i]]) d++; else c++;
        i++;
      }
      if (d > 0) {
        s *= 1 - d / atRisk;
        points.push({ t, s });
      }
      if (c > 0) censors.push({ t, s });
      atRisk -= d + c;
      lastT = t;
    }
    if (points[points.length - 1].t < lastT) points.push({ t: lastT, s });
    return { points, censors };
  }

  function logrank(time1, event1, time2, event2) { // two-group
    const rows = [];
    for (let i = 0; i < time1.length; i++) rows.push([time1[i], event1[i], 0]);
    for (let i = 0; i < time2.length; i++) rows.push([time2[i], event2[i], 1]);
    rows.sort((a, b) => a[0] - b[0]);
    let n = rows.length, n1 = time1.length;
    let O1 = 0, E1 = 0, V = 0;
    let i = 0;
    while (i < rows.length) {
      const t = rows[i][0];
      let d = 0, d1 = 0, rem = 0, rem1 = 0;
      let j = i;
      while (j < rows.length && rows[j][0] === t) {
        if (rows[j][1]) { d++; if (rows[j][2] === 0) d1++; }
        rem++; if (rows[j][2] === 0) rem1++;
        j++;
      }
      if (d > 0) {
        O1 += d1;
        E1 += d * n1 / n;
        if (n > 1) V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1);
      }
      n -= rem; n1 -= rem1;
      i = j;
    }
    const stat = V > 0 ? (O1 - E1) ** 2 / V : 0;
    return { stat, p: chi2Sf(stat, 1) };
  }

  // Cox PH, single covariate, Efron ties (matches lifelines default).
  function coxEfron(time, event, x) {
    // z-score covariate -> HR per 1 SD
    const s = sd(x), m = mean(x);
    if (s === 0) return null;
    const z = x.map((v) => (v - m) / s);
    const idx = time.map((t, i) => i).sort((a, b) => time[a] - time[b]);

    function llGradHess(b) {
      let ll = 0, grad = 0, hess = 0;
      // iterate event times descending; risk set accumulates from the end
      let S0 = 0, S1 = 0, S2 = 0;
      let k = idx.length - 1;
      let i = idx.length - 1;
      while (i >= 0) {
        const t = time[idx[i]];
        while (k >= 0 && time[idx[k]] >= t) {
          const e = Math.exp(b * z[idx[k]]);
          S0 += e; S1 += e * z[idx[k]]; S2 += e * z[idx[k]] ** 2;
          k--;
        }
        let D0 = 0, D1 = 0, D2 = 0, d = 0, sumzx = 0;
        let j = i;
        while (j >= 0 && time[idx[j]] === t) {
          if (event[idx[j]]) {
            const e = Math.exp(b * z[idx[j]]);
            D0 += e; D1 += e * z[idx[j]]; D2 += e * z[idx[j]] ** 2;
            d++; sumzx += z[idx[j]];
          }
          j--;
        }
        if (d > 0) {
          ll += b * sumzx;
          grad += sumzx;
          for (let l = 0; l < d; l++) {
            const f = l / d;
            const den = S0 - f * D0, n1 = S1 - f * D1, n2 = S2 - f * D2;
            ll -= Math.log(den);
            grad -= n1 / den;
            hess -= n2 / den - (n1 / den) ** 2;
          }
        }
        i = j;
      }
      return { ll, grad, hess };
    }

    let b = 0, res = llGradHess(0);
    const ll0 = res.ll;
    for (let it = 0; it < 50; it++) {
      const step = res.grad / -res.hess;
      if (!isFinite(step)) return null;
      b += step;
      res = llGradHess(b);
      if (Math.abs(step) < 1e-10) break;
    }
    const se = 1 / Math.sqrt(-res.hess);
    const lr = 2 * (res.ll - ll0);
    return {
      beta: b, hr: Math.exp(b), se,
      hrLo: Math.exp(b - 1.959963984540054 * se),
      hrHi: Math.exp(b + 1.959963984540054 * se),
      lrStat: lr, lrP: chi2Sf(Math.max(lr, 0), 1),
      waldP: 2 * normSf(Math.abs(b / se)),
      n: time.length, events: event.reduce((a, v) => a + v, 0),
    };
  }

  // maximally selected log-rank cutpoint, 10th-90th pct, >= minPerArm per arm
  function optimalCutpoint(time, event, x, minPerArm = 10) {
    const sorted = [...x].sort((a, b) => a - b);
    const lo = sorted[Math.floor(0.10 * (sorted.length - 1))];
    const hi = sorted[Math.ceil(0.90 * (sorted.length - 1))];
    let candidates = [...new Set(x.filter((v) => v >= lo && v <= hi))].sort((a, b) => a - b);
    if (candidates.length > 150) { // grid-thin
      const th = [];
      for (let i = 0; i < 150; i++) th.push(candidates[Math.round(i * (candidates.length - 1) / 149)]);
      candidates = [...new Set(th)];
    }
    let best = null;
    for (const c of candidates) {
      const hiIdx = [], loIdx = [];
      x.forEach((v, i) => (v > c ? hiIdx : loIdx).push(i));
      if (hiIdx.length < minPerArm || loIdx.length < minPerArm) continue;
      const r = logrank(hiIdx.map((i) => time[i]), hiIdx.map((i) => event[i]),
        loIdx.map((i) => time[i]), loIdx.map((i) => event[i]));
      if (!best || r.stat > best.stat) {
        best = { cutpoint: c, stat: r.stat, logrankP: r.p,
          nHigh: hiIdx.length, nLow: loIdx.length };
      }
    }
    return best;
  }

  globalThis.ChordStats = {
    chi2Sf, tSf2, normSf, mean, median, sd, ranks,
    pearson, spearman, kruskal, kmCurve, logrank, coxEfron, optimalCutpoint,
  };
})();
