## What this changes

<!-- One or two sentences. The diff shows what; say why. -->

## Evidence

<!-- For a behaviour change: which test proves it, and would that test have
     failed before? For a threshold change: on which genomes, with what effect? -->

## Checklist

- [ ] `python -m pytest tests/ -q` passes (full suite, external tools included)
- [ ] Behaviour changes come with a test that would fail without them
- [ ] Any verdict this touches still ships its reasons alongside it
- [ ] Nothing new is random without a fixed seed
- [ ] README, CHANGELOG and the skill file updated if any number or claim moved
