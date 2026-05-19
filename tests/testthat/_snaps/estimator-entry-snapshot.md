# entry-point error: pdata[[time_var]] not integer (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      Invalid inputs:
        - the time_var column 'time' must be integer; got numeric

# entry-point error: time_var not a column (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      Invalid inputs:
        - time_var = 'nonexistent' is not a column in pdata; columns are: unit, time, treat, y

# entry-point error: time_var not character (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      Invalid inputs:
        - time_var must be a single character string; got integer (length 1)

# entry-point error: verbose not logical (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      Invalid inputs:
        - verbose must be a single logical (TRUE or FALSE); got character (length 1)

# entry-point error: alpha out of range (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      Invalid inputs:
        - alpha must be in (0, 1); got -0.1

# entry-point error: fetwfe + betwfe q out of range (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      Invalid inputs:
        - q must be in (0, 2]; got -1

