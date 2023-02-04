mkdir -p ~/.streamlit/

echo "[theme]
base='dark'
primaryColor = ‘#f50c0c’
backgroundColor = ‘#000000’
secondaryBackgroundColor = ‘#bf3c3c’
textColor= ‘#424242’
font = ‘sans serif’
[server]
headless = true
port = $PORT
enableCORS = false
" > ~/.streamlit/config.toml



