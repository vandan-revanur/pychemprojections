from collections import deque


def get_index_of_corresponding_bracket(
    s, i, start_bracket_type="(", end_bracket_type=")"
):
    # If input is invalid.
    if s[i] != start_bracket_type:
        return -1

    # Create a deque to use it as a stack.
    d = deque()

    # Traverse through all elements
    # starting from i.
    for k in range(i, len(s)):
        # Pop a starting bracket
        # for every closing bracket
        if s[k] == end_bracket_type:
            d.popleft()

        # Push all starting brackets
        elif s[k] == start_bracket_type:
            d.append(s[i])

        # If deque becomes empty
        if not d:
            return k

    return -1
