"""Detection codebooks (called languages in the code) built from TE and reference assets."""

from .builder import BUILDER_VERSION, LanguageBuild, build_languages, g1a_tables

__all__ = ["BUILDER_VERSION", "LanguageBuild", "build_languages", "g1a_tables"]
