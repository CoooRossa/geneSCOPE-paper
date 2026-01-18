#!/usr/bin/env bash
set -euo pipefail

OUT=""
SAMPLE_SEC="1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) OUT="$2"; shift 2;;
    --sample_sec) SAMPLE_SEC="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

if [[ -z "${OUT}" ]]; then
  echo "Usage: monitor_cgroup.sh --out <path.tsv> --sample_sec <int>" >&2
  exit 2
fi

SAMPLE_SEC="$(printf "%s" "${SAMPLE_SEC}" | tr -d ' ')"
if [[ -z "${SAMPLE_SEC}" ]]; then SAMPLE_SEC="1"; fi
if [[ "${SAMPLE_SEC}" == "0" ]]; then exit 0; fi

mkdir -p "$(dirname "${OUT}")"

if [[ ! -f "${OUT}" ]]; then
  echo -e "t_sec\tcpu_usage_usec\tmem_bytes" > "${OUT}"
fi

read_cpu_usage() {
  if [[ -r /sys/fs/cgroup/cpu.stat ]]; then
    awk '$1=="usage_usec"{print $2}' /sys/fs/cgroup/cpu.stat
  else
    echo ""
  fi
}

read_mem_bytes() {
  if [[ -r /sys/fs/cgroup/memory.current ]]; then
    cat /sys/fs/cgroup/memory.current
  elif [[ -r /sys/fs/cgroup/memory/memory.usage_in_bytes ]]; then
    cat /sys/fs/cgroup/memory/memory.usage_in_bytes
  else
    echo ""
  fi
}

start_epoch="$(date +%s)"

while true; do
  now_epoch="$(date +%s)"
  t_sec="$((now_epoch - start_epoch))"
  cpu_usec="$(read_cpu_usage || true)"
  mem_bytes="$(read_mem_bytes || true)"
  echo -e "${t_sec}\t${cpu_usec}\t${mem_bytes}" >> "${OUT}"
  sleep "${SAMPLE_SEC}"
done
