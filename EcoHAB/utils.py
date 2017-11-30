def results_path(path):
    head, tail = os.path.split(path)
    head, tail2 = os.path.split(head)
    head = os.path.join(head,'Results_'+tail2)
    return os.path.join(head,tail)
