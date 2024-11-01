

cpsources = custom_target('cpsources',
  command:['cpfiles.sh',meson.current_source_dir(), meson.current_build_dir()],
  output: ['none'],
  build_by_default: true,
)

targets = custom_target('make',
  command:['make', '-C', '@OUTDIR@', 'all'],
  output: ['ddscat', 'calltarget', 'ddpostprocess', 'vtrconvert'],
  build_by_default: true,
  depends: cpsources,
  install: true,
  install_dir: [py.get_install_dir(subdir: 'pyddscat')]
)
