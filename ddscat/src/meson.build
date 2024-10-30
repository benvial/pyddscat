

targets = custom_target('cp',
  command:['cpfiles.sh',meson.current_source_dir(), meson.current_build_dir()],
  output: ['none'],
)

targets = custom_target('make',
  command:['make', '-C', '@OUTDIR@', 'all'],
  output: ['ddscat', 'calltarget', 'ddpostprocess', 'vtrconvert'],
  build_by_default: true,
  depends: targets,
)

