"""Keep validation and research code outside the evitea package."""

import ast
from pathlib import Path


SOURCE_ROOT = Path(__file__).resolve().parents[2] / "src" / "evitea"


def _is_forbidden_module(module: str) -> bool:
    parts = [part.casefold() for part in module.split(".")]
    joined = ".".join(parts)
    return (
        "validation" in parts
        or any("implementation_prep" in part for part in parts)
        or ("benchmark" in joined and "truth" in joined)
        or any("temp2" in part or "teforest" in part for part in parts)
    )


def _imported_modules(tree: ast.AST) -> list[tuple[str, int]]:
    imports: list[tuple[str, int]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend((alias.name, node.lineno) for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                imports.append((node.module, node.lineno))
            imports.extend((alias.name, node.lineno) for alias in node.names)
        elif isinstance(node, ast.Call):
            function = node.func
            is_dynamic_import = (
                isinstance(function, ast.Name) and function.id == "__import__"
            ) or (
                isinstance(function, ast.Attribute)
                and function.attr == "import_module"
                and isinstance(function.value, ast.Name)
                and function.value.id == "importlib"
            )
            if is_dynamic_import and node.args:
                argument = node.args[0]
                if isinstance(argument, ast.Constant) and isinstance(argument.value, str):
                    imports.append((argument.value, node.lineno))
    return imports


def test_production_tree_has_no_validation_or_research_history_imports() -> None:
    violations: list[str] = []
    for path in sorted(SOURCE_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for module, line in _imported_modules(tree):
            if _is_forbidden_module(module):
                violations.append(f"{path.relative_to(SOURCE_ROOT)}:{line}: {module}")

    assert not violations, "Forbidden production imports:\n" + "\n".join(violations)
