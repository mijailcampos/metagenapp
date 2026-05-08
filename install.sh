#!/bin/bash
set -e

CYAN='\033[1;36m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
RESET='\033[0m'

echo ""
echo -e "${BOLD}${CYAN}███╗   ███╗███████╗████████╗ █████╗  ██████╗ ███████╗███╗   ██╗${RESET}"
echo -e "${BOLD}${CYAN}████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝████╗  ██║${RESET}"
echo -e "${BOLD}${CYAN}██╔████╔██║█████╗     ██║   ███████║██║  ███╗█████╗  ██╔██╗ ██║${RESET}"
echo -e "${BOLD}${CYAN}██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██╔══╝  ██║╚██╗██║${RESET}"
echo -e "${BOLD}${CYAN}██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝███████╗██║ ╚████║${RESET}"
echo -e "${BOLD}${CYAN}╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═══╝${RESET}"
echo ""
echo -e "           ${BOLD}Installer — Metabarcoding 16S/18S${RESET}"
echo ""

# ── 1. Verificar conda ──────────────────────────────────────────
echo -e "${CYAN}[1/4] Verificando conda...${RESET}"
if ! command -v conda &>/dev/null; then
    echo -e "${RED}✗ conda no encontrado.${RESET}"
    echo -e "  Instala Miniconda desde: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi
echo -e "${GREEN}  ✓ conda $(conda --version)${RESET}"

# ── 2. Crear entorno ────────────────────────────────────────────
echo ""
echo -e "${CYAN}[2/4] Creando entorno conda (metagenapp)...${RESET}"
echo -e "  ${YELLOW}Esto puede tardar unos minutos.${RESET}"
conda env create -f environment.yml

# ── 3. Descargar modelo SILVA v1 ────────────────────────────────
echo ""
echo -e "${CYAN}[3/4] Descargando modelo SILVA v1 (~2.3 GB)...${RESET}"
echo -e "  ${YELLOW}Esto puede tardar varios minutos dependiendo de tu conexión.${RESET}"
conda run -n metagenapp python download_model.py

# ── 4. Listo ────────────────────────────────────────────────────
echo ""
echo -e "${GREEN}${BOLD}✓ Instalación completa.${RESET}"
echo ""
echo -e "  Para usar MetagenApp:"
echo -e ""
echo -e "    ${BOLD}conda activate metagenapp${RESET}"
echo -e "    ${BOLD}metagenapp-launch${RESET}       ← menú interactivo"
echo -e "    ${BOLD}metagenapp --help${RESET}        ← opciones del CLI"
echo ""
