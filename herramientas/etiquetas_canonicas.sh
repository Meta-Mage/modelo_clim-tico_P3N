#!/usr/bin/env bash
# herramientas/etiquetas_canonicas.sh -- M3N 3.12.1 (NOMENCLATURA.md §3.1 y §4.1).
# Crea las etiquetas canonicas (SemVer) como ALIAS anotados del MISMO commit que las etiquetas antiguas.
# No borra ni mueve ninguna etiqueta. Es idempotente: si la canonica ya existe en el mismo commit, no hace nada;
# si existe en OTRO commit, para con un error sin tocar nada mas.
# Uso (en la carpeta de M3N):  bash herramientas/etiquetas_canonicas.sh  &&  git push origin --tags
set -euo pipefail
PARES="v2.0:v2.0.0 v2.1:v2.1.0 v2.2:v2.2.0 v2.2b:v2.2.1 v2.3:v2.3.0 v3.0-pre1:v3.0.0
v3.1-pre3:v3.1.0 v3.1-pre4:v3.2.0 v3.1-pre5:v3.3.0 v3.1-pre6:v3.4.0 v3.1-pre7:v3.5.0 v3.1-pre8:v3.6.0
v3.1-pre9:v3.7.0 v3.1-pre10:v3.8.0 v3.1-pre11:v3.9.0 v3.1-pre12:v3.10.0 v3.1-pre13:v3.10.1
v3.1-pre14:v3.10.2 v3.1-pre15:v3.11.0 v3.1-pre16:v3.12.0"
creadas=0; existian=0; sin_antigua=0
for par in $PARES; do
  vieja="${par%%:*}"; nueva="${par##*:}"
  if ! git rev-parse -q --verify "refs/tags/$vieja" >/dev/null; then
    echo "  (no existe la etiqueta antigua $vieja: sin alias)"; sin_antigua=$((sin_antigua + 1)); continue
  fi
  commit="$(git rev-list -n 1 "$vieja")"
  if git rev-parse -q --verify "refs/tags/$nueva" >/dev/null; then
    if [ "$(git rev-list -n 1 "$nueva")" != "$commit" ]; then
      echo "ERROR: $nueva ya existe y apunta a otro commit que $vieja. No sigo."; exit 1
    fi
    existian=$((existian + 1)); continue
  fi
  git tag -a "$nueva" "$commit" -m "M3N ${nueva#v} (alias de $vieja; nomenclatura normalizada, NOMENCLATURA.md)"
  echo "  creada $nueva -> $vieja (${commit:0:7})"; creadas=$((creadas + 1))
done
echo "Etiquetas canonicas: $creadas creadas, $existian ya existian, $sin_antigua antiguas inexistentes."
