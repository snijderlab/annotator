const { openPath } = window.__TAURI__.opener;
const { resolveResource, resourceDir } = window.__TAURI__.path;

window.addEventListener("DOMContentLoaded", () => {
    document.querySelectorAll("a[target=_blank]").forEach(element => element.addEventListener("click", async (e) => {
        e.preventDefault();
        if (e.target.dataset.href != undefined) {
            openPath(await resolveResource(e.target.dataset.href))
        } else {
            openPath(e.target.href);
        }
    }));
})