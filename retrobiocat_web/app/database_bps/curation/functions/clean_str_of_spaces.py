

def remove_end_spaces_if_str(var):
    if not isinstance(var, str):
        return var

    if len(var) == 0:
        return var

    while var[-1] == ' ':
        var = var[:-1]
        if len(var) == 0:
            return var
    return var