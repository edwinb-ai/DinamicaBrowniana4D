project('dbrownian', 'c')
add_project_arguments('-march=native', language : 'c')
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)
executable('db', 
    'langevin_4D.c', 'functionsl.c', 'splitmix64.c', 'xoshiro256plus.c',
    'random_s.c',
    dependencies : m_dep)
