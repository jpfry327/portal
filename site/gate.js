/* Pre-publication access gate.
 *
 * SCOPE — read this before relying on it. This is a deterrent against
 * ACCIDENTAL discovery, not a security control. The site is static and served
 * from a public repository, so anyone willing to open developer tools or read
 * the repo can reach the data regardless of this gate. It exists so that a
 * search engine, a shared link, or a guessed URL does not put unpublished
 * results in front of someone who was not looking for them.
 *
 * For genuine access control, put the site behind Cloudflare Access (free tier
 * covers 50 users) so unauthenticated requests never reach the content.
 *
 * The passphrase is compared as a SHA-256 digest so it is not sitting in the
 * source as plain text — this raises the effort slightly; it does not make the
 * gate strong, since the check runs client-side either way.
 */
(function () {
  // sha256("chordoma")
  const EXPECT = "8881231063d0889f60a113f27292025a8e6c9b2ad2752ccc1084e426c8d8b3a7";
  const KEY = "chordoma-atlas-unlocked";

  async function sha256(s) {
    const buf = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(s));
    return [...new Uint8Array(buf)].map((b) => b.toString(16).padStart(2, "0")).join("");
  }

  function screen() {
    const d = document.createElement("div");
    d.id = "gate";
    d.innerHTML = `
      <div class="gate-card">
        <h2>Access restricted</h2>
        <p>This resource accompanies a manuscript under review and is not yet public.
           Please enter the passphrase to continue.</p>
        <form id="gate-form" autocomplete="off">
          <input type="password" id="gate-input" placeholder="Passphrase"
                 autocomplete="current-password" spellcheck="false" autofocus>
          <button type="submit">Enter</button>
        </form>
        <div id="gate-msg"></div>
      </div>`;
    return d;
  }

  const STYLE = `
    #gate { position:fixed; inset:0; z-index:9999; display:flex; align-items:center;
            justify-content:center; background:linear-gradient(135deg,#440154,#31688E 55%,#35B779);
            font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; }
    .gate-card { background:#fff; border-radius:12px; padding:30px 34px; max-width:420px;
                 box-shadow:0 12px 40px rgba(0,0,0,.3); }
    .gate-card h2 { margin:0 0 10px; font-size:20px; color:#1a1a2e; }
    .gate-card p { margin:0 0 18px; font-size:14px; line-height:1.5; color:#556; }
    #gate-form { display:flex; gap:8px; }
    #gate-input { flex:1; padding:10px 12px; border:1.5px solid #d5d5df; border-radius:7px; font-size:14px; }
    #gate-input:focus { outline:none; border-color:#31688E; }
    #gate-form button { padding:10px 20px; border:none; border-radius:7px; background:#440154;
                        color:#fff; font-size:14px; font-weight:600; cursor:pointer; }
    #gate-msg { margin-top:10px; font-size:13px; color:#b2182b; min-height:18px; }
    body.locked > *:not(#gate) { display:none !important; }`;

  // app.js is a deferred module, so it may register its listener AFTER this
  // file runs. Latch the state on globalThis as well as firing the event, and
  // let app.js check the latch first — otherwise an already-unlocked session
  // fires the event into the void and the app never boots.
  function signalUnlock() {
    globalThis.ATLAS_UNLOCKED = true;
    document.dispatchEvent(new Event("atlas-unlock"));
  }

  function unlock(el) {
    document.body.classList.remove("locked");
    el.remove();
    signalUnlock();
  }

  document.addEventListener("DOMContentLoaded", async () => {
    const style = document.createElement("style");
    style.textContent = STYLE;
    document.head.appendChild(style);

    if (sessionStorage.getItem(KEY) === "1") {
      signalUnlock();
      return;
    }

    document.body.classList.add("locked");
    const el = screen();
    document.body.appendChild(el);
    document.getElementById("gate-input").focus();

    document.getElementById("gate-form").addEventListener("submit", async (e) => {
      e.preventDefault();
      const val = document.getElementById("gate-input").value;
      if (await sha256(val) === EXPECT) {
        sessionStorage.setItem(KEY, "1");
        unlock(el);
      } else {
        document.getElementById("gate-msg").textContent = "Incorrect passphrase.";
        document.getElementById("gate-input").select();
      }
    });
  });
})();
