
def is_overlapping(long x1,long x2,long y1,long y2):
    '''
    assumes x1/y1 are mins among group, returns bool True/False whether there's overlap or not
    '''
    if x1 > x2:
        x_min = x2
        x_max = x1
    else:
        x_min = x1
        x_max = x2

    if y1 > y2:
        y_min = y2
        y_max = y1
    else:
        y_min = y1
        y_max = y2

    if x_min > y_min:
        left_side = x_min
    else:
        left_side = y_min
        
    if x_max > y_max:
        right_side = y_max
    else:
        right_side = x_max

    return left_side <= right_side
#    return max(x_min, y_min) <= min(x_max, y_max)
