from collections import deque


def get_index_of_corresponding_bracket(
    input_string: str,
    start_bracket_idx: int,
    start_bracket_type: str = "(",
    end_bracket_type: str = ")",
) -> int:
    # If input is invalid.
    if input_string[start_bracket_idx] != start_bracket_type:
        return -1

    # Create a deque to use it as a stack.
    d = deque()

    # Traverse through all elements
    # starting from start_bracket_idx.
    for k in range(start_bracket_idx, len(input_string)):
        # Pop a starting bracket
        # for every closing bracket
        if input_string[k] == end_bracket_type:
            d.popleft()

        # Push all starting brackets
        elif input_string[k] == start_bracket_type:
            d.append(input_string[start_bracket_idx])

        # If deque becomes empty
        if not d:
            return k

    return -1
