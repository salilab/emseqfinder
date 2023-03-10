import saliweb.build

vars = Variables('config.py')
env = saliweb.build.Environment(vars, ['conf/live.conf'], service_module='emseqfinder')
Help(vars.GenerateHelpText(env))

env.InstallAdminTools()

Export('env')
SConscript('frontend/emseqfinder/SConscript')
SConscript('backend/emseqfinder/SConscript')
SConscript('test/SConscript')
