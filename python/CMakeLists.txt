file(GLOB to_install *py figure.txt)

install(FILES ${to_install} DESTINATION lib/python2.7/site-packages/podio)
