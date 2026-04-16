#!/usr/bin/env bash
# start.sh — Levanta backend (uvicorn) y frontend (vite) de MetagenApp
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cleanup() {
    echo ""
    echo "Deteniendo servicios..."
    kill "$BACKEND_PID" "$FRONTEND_PID" 2>/dev/null
    wait "$BACKEND_PID" "$FRONTEND_PID" 2>/dev/null
    echo "Listo."
}
trap cleanup EXIT INT TERM

echo "Iniciando backend (uvicorn) en :8000 ..."
cd "$SCRIPT_DIR"
uvicorn metagenapp.api.main:app --reload --port 8000 --host 0.0.0.0 &
BACKEND_PID=$!

echo "Iniciando frontend (vite) ..."
cd "$SCRIPT_DIR/frontend"
npm run dev -- --host &
FRONTEND_PID=$!

echo ""
echo "  Backend:  http://0.0.0.0:8000"
echo "  Frontend: http://0.0.0.0:5173  (o el puerto que muestre vite arriba)"
echo ""
echo "Ctrl+C para detener ambos."

wait
