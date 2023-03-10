
async function select_file(e) {
  let properties = {
    //defaultPath: 'C:\\',
    directory: false,
    filters: [{
      extensions: ['mgf'], name: "*"
    }]
  };
  e.dataset.filepath = await window.__TAURI__.dialog.open(properties);
};