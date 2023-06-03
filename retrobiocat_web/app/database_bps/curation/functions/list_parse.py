

def list_to_string(list_field):
    list_string = ""
    for list_item in list_field:
        if list_string != "":
            list_string += ", "
        list_string += list_item
    return list_string