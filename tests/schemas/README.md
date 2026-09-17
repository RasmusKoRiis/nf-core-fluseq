# Vendored module metadata schema

`module-meta.schema.json` is the nf-core/modules [module metadata schema](https://github.com/nf-core/modules/blob/master/modules/meta-schema.json), downloaded on 2026-09-11. The upstream MIT license is included as `nf-core-modules-LICENSE.txt`.

The SHA-256 of the original downloaded JSON, before repository formatting, was `9458d6897d70abb4aa20d207cdf7ddca5162ffd48f3b71827ec409d51f087051`.

The copy keeps metadata tests independent of network access and upstream changes. To update it, review the upstream schema changes, replace the copy, retain attribution, and run `pytest tests/test_module_metadata.py`. Schema validation checks documentation structure; it does not validate biological output or software environments.
