# test_nomenclatura.py -- M3N 3.12.1: el repositorio cumple la norma de NOMENCLATURA.md (versiones SemVer, CHANGELOG,
# documentos de diseño por fase, sin nombres antiguos fuera de las tablas de correspondencia, interruptores).
#
# Uso:  python -m pytest test_nomenclatura.py      (segundos)

import os
import re

AQUI = os.path.dirname(os.path.abspath(__file__))
SEMVER = re.compile(r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[A-Za-z-][0-9A-Za-z-]*)"
                    r"(?:\.(?:0|[1-9]\d*|\d*[A-Za-z-][0-9A-Za-z-]*))*))?$")       # semver.org, sin metadatos
# archivos donde SÍ pueden aparecer los nombres antiguos (son las tablas de correspondencia)
CORRESPONDENCIAS = {"NOMENCLATURA.md", "CHANGELOG.md", os.path.join("herramientas", "etiquetas_canonicas.sh"),
                    "test_nomenclatura.py"}
DOCS_ANTIGUOS = ("NOTAS_DISENO_FASE0.md", "DISENO_FASE2B.md", "DISENO_V3.0.md", "DISENO_FASE3.md", "DISENO_FASE4.md",
                 "DISENO_FASE5A.md", "DISENO_V3.1.md", "DISENO_FASE6_BORRADOR.md", "DISENO_FASE6_1.md",
                 "DISENO_FASE6_2.md", "DISENO_FASE6_3.md")


def _textos():
    """Archivos de texto del repositorio (sin outputs, venv ni cachés)."""
    for raiz, dirs, archivos in os.walk(AQUI):
        dirs[:] = [d for d in dirs if d not in ("outputs", "venv", "__pycache__", ".git", ".pytest_cache", "bench")]
        for a in archivos:
            if a.endswith((".py", ".md", ".sh", ".txt")):
                ruta = os.path.join(raiz, a)
                yield os.path.relpath(ruta, AQUI), open(ruta, encoding="utf-8").read()


def _lineas_vigentes(rel, texto):
    """Líneas sujetas a la norma: fuera de las tablas de correspondencia y de las notas de nomenclatura."""
    if rel in CORRESPONDENCIAS:
        return []
    return [(i + 1, l) for i, l in enumerate(texto.split("\n")) if not l.startswith("> **Nomenclatura")]


def _semver_clave(v):
    m = SEMVER.match(v)
    nucleo = tuple(int(x) for x in m.group(1, 2, 3))
    pre = m.group(4)
    if pre is None:
        return nucleo + ((1,),)                                  # una versión final va después de sus previas
    ids = tuple((0, int(x), "") if x.isdigit() else (1, 0, x) for x in pre.split("."))
    return nucleo + ((0,) + ids,)


def test_version_unica_y_semver():
    """VERSION (clima_dinamico.py), la línea del README y la primera entrada del CHANGELOG coinciden y son SemVer."""
    codigo = open(os.path.join(AQUI, "clima_dinamico.py"), encoding="utf-8").read()
    v_codigo = re.search(r'^VERSION = "([^"]+)"', codigo, re.M).group(1)
    v_readme = re.search(r"Versión actual: `([^`]+)`", open(os.path.join(AQUI, "README.md"), encoding="utf-8").read()).group(1)
    cambios = open(os.path.join(AQUI, "CHANGELOG.md"), encoding="utf-8").read()
    v_changelog = re.search(r"^## \[([^\]]+)\]", cambios, re.M).group(1)
    assert SEMVER.match(v_codigo), v_codigo
    assert v_codigo == v_readme == v_changelog, (v_codigo, v_readme, v_changelog)


def test_changelog_ordenado_y_valido():
    """Todas las versiones del CHANGELOG son SemVer, sin repetir, de la más reciente a la más antigua, y cada una
    dice su fase."""
    cambios = open(os.path.join(AQUI, "CHANGELOG.md"), encoding="utf-8").read()
    versiones = re.findall(r"^## \[([^\]]+)\]", cambios, re.M)
    assert len(versiones) >= 25
    for v in versiones:
        assert SEMVER.match(v), v
    assert len(set(versiones)) == len(versiones)
    claves = [_semver_clave(v) for v in versiones]
    assert claves == sorted(claves, reverse=True), "el CHANGELOG no va de la más reciente a la más antigua"
    bloques = re.split(r"^## \[", cambios, flags=re.M)[1:]
    for b in bloques:
        assert "\n- Fase:" in b, b.split("\n")[0]


def test_semver_ordena_bien_y_pre_no():
    """La razón de la norma: con SemVer, 'pre16' iría antes que 'pre6'; con identificadores separados por punto, no."""
    assert _semver_clave("3.1.0-pre16") < _semver_clave("3.1.0-pre6")         # el defecto de los nombres antiguos
    assert _semver_clave("4.0.0-rc.2") < _semver_clave("4.0.0-rc.10") < _semver_clave("4.0.0")
    assert _semver_clave("3.10.0") > _semver_clave("3.9.0") > _semver_clave("3.1.0")


def test_documentos_de_diseno():
    """Cada documento de diseño se llama DISENO_FASE_<id>.md con un id de fase válido, está en la tabla del README
    y en la de NOMENCLATURA; cada DISENO_FASE_*.md citado existe; los nombres antiguos ya no existen."""
    docs = sorted(f for f in os.listdir(AQUI) if f.startswith(("DISENO", "NOTAS_DISENO")) and f.endswith(".md"))
    assert docs, "no hay documentos de diseño"
    readme = open(os.path.join(AQUI, "README.md"), encoding="utf-8").read()
    norma = open(os.path.join(AQUI, "NOMENCLATURA.md"), encoding="utf-8").read()
    for d in docs:
        assert re.fullmatch(r"DISENO_FASE_\d+(\.\d+)*\.md", d), d
        assert f"`{d}`" in readme, f"{d} no está en la tabla del README"
        assert f"`{d}`" in norma, f"{d} no está en NOMENCLATURA.md"
    for viejo in DOCS_ANTIGUOS:
        assert not os.path.exists(os.path.join(AQUI, viejo)), viejo
    for rel, texto in _textos():
        for citado in set(re.findall(r"DISENO_FASE_[\d.]*\d\.md", texto)):
            assert os.path.exists(os.path.join(AQUI, citado)), f"{rel} cita {citado}, que no existe"


def test_sin_nombres_antiguos():
    """Fuera de las tablas de correspondencia no quedan nombres antiguos de documentos, versiones ni fases."""
    patrones = [re.compile(re.escape(d)) for d in DOCS_ANTIGUOS] + [
        re.compile(r"v3\.[01]-pre\d"), re.compile(r"(?<![\w.-])pre\d+\b"), re.compile(r"\bv2\.2[bc]\b"),
        re.compile(r"\b[Ff]ases? \d+[a-z]\b"), re.compile(r"\bla \d[ab]\b")]
    malos = []
    for rel, texto in _textos():
        for n, linea in _lineas_vigentes(rel, texto):
            for p in patrones:
                if p.search(linea):
                    malos.append(f"{rel}:{n}: {p.pattern}: {linea.strip()[:100]}")
    assert not malos, "\n".join(malos[:40])


def test_interruptores_documentados():
    """La tabla de interruptores de NOMENCLATURA §6 tiene exactamente los del código, con su número y en orden."""
    import fase2b_atmosfera as F
    norma = open(os.path.join(AQUI, "NOMENCLATURA.md"), encoding="utf-8").read()
    filas = re.findall(r"^\| I(\d+) \| `([a-z_]+)` \|", norma, re.M)
    assert [int(n) for n, _ in filas] == list(range(1, len(F.INTERRUPTORES_FASE2B) + 1))
    assert [k for _, k in filas] == list(F.INTERRUPTORES_FASE2B)


def test_herramientas():
    """Las herramientas auxiliares viven en herramientas/ y el script de etiquetas cubre todas las antiguas."""
    sh = open(os.path.join(AQUI, "herramientas", "etiquetas_canonicas.sh"), encoding="utf-8").read()
    pares = dict(p.split(":") for p in re.findall(r"v[\d.]+(?:-pre\d+|b)?:v\d+\.\d+\.\d+", sh))
    for n in range(3, 17):
        assert f"v3.1-pre{n}" in pares
    assert all(SEMVER.match(v[1:]) for v in pares.values())
    assert len(set(pares.values())) == len(pares)
    assert os.path.exists(os.path.join(AQUI, "herramientas", "diagnostico_tierra.py"))
