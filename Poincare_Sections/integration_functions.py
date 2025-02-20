base_funcs = {
    "b$110$000$000": 'Simplify[D[Integrate[1, {x, left, l2}, {y, func, bottom1}], t]]',
    "b$010$000$000": 'Simplify[D[Integrate[1, {x, l1, l2}, {y, func, bottom1}], t]]',
    "b$011$000$000": 'Simplify[D[Integrate[1, {x, l1, point1}, {y, func, bottom1}], t]]',
    "b$111$000$000": 'Simplify[D[Integrate[1, {x, left, point1}, {y, func, bottom1}], t]]',

    "b$000$110$000": 'Simplify[D[Integrate[1, {x, point1, m2}, {y, func, bottom2}], t]]',
    "b$000$010$000": 'Simplify[D[Integrate[1, {x, m1, m2}, {y, func, bottom2}], t]]',
    "b$000$011$000": 'Simplify[D[Integrate[1, {x, m1, point2}, {y, func, bottom2}], t]]',
    "b$000$111$000": 'Simplify[D[Integrate[1, {x, point1, point2}, {y, func, bottom2}], t]]',

    "b$000$000$110": 'Simplify[D[Integrate[1, {x, point2, r2}, {y, func, bottom3}], t]]',
    "b$000$000$010": 'Simplify[D[Integrate[1, {x, r1, r2}, {y, func, bottom3}], t]]',
    "b$000$000$011": 'Simplify[D[Integrate[1, {x, r1, 1}, {y, func, bottom3}], t]]',
    "b$000$000$111": 'Simplify[D[Integrate[1, {x, point2, 1}, {y, func, bottom3}], t]]'
}

equations_dict = {
    
#one edge

    # bottom1
    "f$110$000$000": f"f$110$000$000 = Simplify[fL - {base_funcs['b$110$000$000']}]",
    "f$010$000$000": f"f$010$000$000 = Simplify[f1 - {base_funcs['b$010$000$000']}]",
    "f$011$000$000": f"f$011$000$000 = Simplify[f1 - {base_funcs['b$011$000$000']}]",
    "f$111$000$000": f"f$111$000$000 = Simplify[fL - {base_funcs['b$111$000$000']}]",

    # bottom2
    "f$000$110$000": f"f$000$110$000 = Simplify[f1 - {base_funcs['b$000$110$000']}]",
    "f$000$010$000": f"f$000$010$000 = Simplify[f1 - {base_funcs['b$000$010$000']}]",
    "f$000$011$000": f"f$000$011$000 = Simplify[f1 - {base_funcs['b$000$011$000']}]",
    "f$000$111$000": f"f$000$111$000 = Simplify[f1 - {base_funcs['b$000$111$000']}]",

    # bottom3
    "f$000$000$110": f"f$000$000$110 = Simplify[f1 - {base_funcs['b$000$000$110']}]",
    "f$000$000$010": f"f$000$000$010 = Simplify[f1 - {base_funcs['b$000$000$010']}]",
    "f$000$000$011": f"f$000$000$011 = Simplify[f1 - {base_funcs['b$000$000$011']}]",
    "f$000$000$111": f"f$000$000$111 = Simplify[f1 - {base_funcs['b$000$000$111']}]",

#################################################################################################################################

# two edges

    # bottom1 110
        # bottom2
    "f$110$110$000": f"f$110$110$000 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$110$000']}]",
    "f$110$010$000": f"f$110$010$000 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$010$000']}]",
    "f$110$011$000": f"f$110$011$000 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$011$000']}]",
    "f$110$111$000": f"f$110$111$000 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$111$000']}]",
        # bottom3
    "f$110$000$110": f"f$110$000$110 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$000$110']}]",
    "f$110$000$010": f"f$110$000$010 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$000$010']}]",
    "f$110$000$011": f"f$110$000$011 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$000$011']}]",
    "f$110$000$111": f"f$110$000$111 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$000$111']}]",

    # bottom1 010
        # bottom2
    "f$010$110$000": f"f$010$110$000 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$110$000']}]",
    "f$010$010$000": f"f$010$010$000 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$010$000']}]",
    "f$010$011$000": f"f$010$011$000 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$011$000']}]",
    "f$010$111$000": f"f$010$111$000 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$111$000']}]",
        # bottom3
    "f$010$000$110": f"f$010$000$110 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$000$110']}]",
    "f$010$000$010": f"f$010$000$010 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$000$010']}]",
    "f$010$000$011": f"f$010$000$011 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$000$011']}]",
    "f$010$000$111": f"f$010$000$111 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$000$111']}]",

    # bottom1 011
        # bottom2
    "f$011$110$000": f"f$011$110$000 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$110$000']}]",
    "f$011$010$000": f"f$011$010$000 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$010$000']}]",
    "f$011$011$000": f"f$011$011$000 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$011$000']}]",
    "f$011$111$000": f"f$011$111$000 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$111$000']}]",
        # bottom3
    "f$011$000$110": f"f$011$000$110 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$000$110']}]",
    "f$011$000$010": f"f$011$000$010 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$000$010']}]",
    "f$011$000$011": f"f$011$000$011 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$000$011']}]",
    "f$011$000$111": f"f$011$000$111 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$000$111']}]",

    # bottom1 111
        # bottom2
    "f$111$110$000": f"f$111$110$000 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$110$000']}]",
    "f$111$010$000": f"f$111$010$000 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$010$000']}]",
    "f$111$011$000": f"f$111$011$000 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$011$000']}]",
    "f$111$111$000": f"f$111$111$000 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$111$000']}]",
        # bottom3
    "f$111$000$110": f"f$111$000$110 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$000$110']}]",
    "f$111$000$010": f"f$111$000$010 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$000$010']}]",
    "f$111$000$011": f"f$111$000$011 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$000$011']}]",
    "f$111$000$111": f"f$111$000$111 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$000$111']}]",

    # bottom2 110
        # bottom3
    "f$000$110$110": f"f$000$110$110 = Simplify[f1 - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$110']}]",
    "f$000$110$010": f"f$000$110$010 = Simplify[f1 - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$010']}]",
    "f$000$110$011": f"f$000$110$011 = Simplify[f1 - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$011']}]",
    "f$000$110$111": f"f$000$110$111 = Simplify[f1 - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$111']}]",
   
    # bottom2 010
        # bottom3
    "f$000$010$110": f"f$000$010$110 = Simplify[f1 - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$110']}]",
    "f$000$010$010": f"f$000$010$010 = Simplify[f1 - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$010']}]",
    "f$000$010$011": f"f$000$010$011 = Simplify[f1 - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$011']}]",
    "f$000$010$111": f"f$000$010$111 = Simplify[f1 - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$111']}]",

    # bottom2 011
        # bottom3
    "f$000$011$110": f"f$000$011$110 = Simplify[f1 - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$110']}]",
    "f$000$011$010": f"f$000$011$010 = Simplify[f1 - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$010']}]",
    "f$000$011$011": f"f$000$011$011 = Simplify[f1 - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$011']}]",
    "f$000$011$111": f"f$000$011$111 = Simplify[f1 - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$111']}]",

    # bottom2 111
        # bottom3
    "f$000$111$110": f"f$000$111$110 = Simplify[f1 - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$110']}]",
    "f$000$111$010": f"f$000$111$010 = Simplify[f1 - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$010']}]",
    "f$000$111$011": f"f$000$111$011 = Simplify[f1 - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$011']}]",
    "f$000$111$111": f"f$000$111$111 = Simplify[f1 - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$111']}]",

#################################################################################################################################

# three edges

    # bottom1 110
        # bottom2 110
            # bottom3
    "f$110$110$110": f"f$110$110$110 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$110']}]",
    "f$110$110$010": f"f$110$110$010 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$010']}]",
    "f$110$110$011": f"f$110$110$011 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$011']}]",
    "f$110$110$111": f"f$110$110$111 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 010
            # bottom3
    "f$110$010$110": f"f$110$010$110 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$110']}]",
    "f$110$010$010": f"f$110$010$010 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$010']}]",
    "f$110$010$011": f"f$110$010$011 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$011']}]",
    "f$110$010$111": f"f$110$010$111 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 011
            # bottom3
    "f$110$011$110": f"f$110$011$110 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$110']}]",
    "f$110$011$010": f"f$110$011$010 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$010']}]",
    "f$110$011$011": f"f$110$011$011 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$011']}]",
    "f$110$011$111": f"f$110$011$111 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 111
            # bottom3
    "f$110$111$110": f"f$110$111$110 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$110']}]",
    "f$110$111$010": f"f$110$111$010 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$010']}]",
    "f$110$111$011": f"f$110$111$011 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$011']}]",
    "f$110$111$111": f"f$110$111$111 = Simplify[fL - {base_funcs['b$110$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$111']}]",

    #bottom1 010
        # bottom2 110
            # bottom3
    "f$010$110$110": f"f$010$110$110 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$110']}]",
    "f$010$110$010": f"f$010$110$010 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$010']}]",
    "f$010$110$011": f"f$010$110$011 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$011']}]",
    "f$010$110$111": f"f$010$110$111 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 010
            # bottom3
    "f$010$010$110": f"f$010$010$110 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$110']}]",
    "f$010$010$010": f"f$010$010$010 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$010']}]",
    "f$010$010$011": f"f$010$010$011 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$011']}]",
    "f$010$010$111": f"f$010$010$111 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 011
            # bottom3
    "f$010$011$110": f"f$010$011$110 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$110']}]",
    "f$010$011$010": f"f$010$011$010 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$010']}]",
    "f$010$011$011": f"f$010$011$011 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$011']}]",
    "f$010$011$111": f"f$010$011$111 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 111
            # bottom3
    "f$010$111$110": f"f$010$111$110 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$110']}]",
    "f$010$111$010": f"f$010$111$010 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$010']}]",
    "f$010$111$011": f"f$010$111$011 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$011']}]",
    "f$010$111$111": f"f$010$111$111 = Simplify[f1 - {base_funcs['b$010$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$111']}]",

    #bottom1 011
        # bottom2 110
            # bottom3
    "f$011$110$110": f"f$011$110$110 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$110']}]",
    "f$011$110$010": f"f$011$110$010 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$010']}]",
    "f$011$110$011": f"f$011$110$011 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$011']}]",
    "f$011$110$111": f"f$011$110$111 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 010
            # bottom3
    "f$011$010$110": f"f$011$010$110 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$110']}]",
    "f$011$010$010": f"f$011$010$010 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$010']}]",
    "f$011$010$011": f"f$011$010$011 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$011']}]",
    "f$011$010$111": f"f$011$010$111 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 011
            # bottom3
    "f$011$011$110": f"f$011$011$110 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$110']}]",
    "f$011$011$010": f"f$011$011$010 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$010']}]",
    "f$011$011$011": f"f$011$011$011 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$011']}]",
    "f$011$011$111": f"f$011$011$111 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 111
            # bottom3
    "f$011$111$110": f"f$011$111$110 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$110']}]",
    "f$011$111$010": f"f$011$111$010 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$010']}]",
    "f$011$111$011": f"f$011$111$011 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$011']}]",
    "f$011$111$111": f"f$011$111$111 = Simplify[f1 - {base_funcs['b$011$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$111']}]",

    #bottom1 111
        # bottom2 110
            # bottom3
    "f$111$110$110": f"f$111$110$110 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$110']}]",
    "f$111$110$010": f"f$111$110$010 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$010']}]",
    "f$111$110$011": f"f$111$110$011 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$011']}]",
    "f$111$110$111": f"f$111$110$111 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$110$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 010
            # bottom3
    "f$111$010$110": f"f$111$010$110 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$110']}]",
    "f$111$010$010": f"f$111$010$010 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$010']}]",
    "f$111$010$011": f"f$111$010$011 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$011']}]",
    "f$111$010$111": f"f$111$010$111 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$010$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 011
            # bottom3
    "f$111$011$110": f"f$111$011$110 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$110']}]",
    "f$111$011$010": f"f$111$011$010 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$010']}]",
    "f$111$011$011": f"f$111$011$011 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$011']}]",
    "f$111$011$111": f"f$111$011$111 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$011$000']} - {base_funcs['b$000$000$111']}]",

        # bottom2 111
            # bottom3
    "f$111$111$110": f"f$111$111$110 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$110']}]",
    "f$111$111$010": f"f$111$111$010 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$010']}]",
    "f$111$111$011": f"f$111$111$011 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$011']}]",
    "f$111$111$111": f"f$111$111$111 = Simplify[fL - {base_funcs['b$111$000$000']} - {base_funcs['b$000$111$000']} - {base_funcs['b$000$000$111']}]"
}