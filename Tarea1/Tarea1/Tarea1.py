import re
cadena = 'GATTATATACATAGTAGTATA'

exp = re.compile('(TA)+')

for match in re.finditer(exp, cadena):
    s = match.start()
    f = match.end()
    print('Se encontró %s en la posicion %d:%d' % (cadena[s:f], s, f))