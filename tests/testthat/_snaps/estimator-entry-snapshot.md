# entry-point error: pdata[[time_var]] not integer (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      is.integer(pdata[[time_var]]) is not TRUE

# entry-point error: time_var not a column (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      time_var %in% colnames(pdata) is not TRUE

# entry-point error: time_var not character (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      is.character(time_var) is not TRUE

# entry-point error: verbose not logical (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      is.logical(verbose) is not TRUE

# entry-point error: alpha out of range (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      alpha > 0 is not TRUE

# entry-point error: fetwfe + betwfe q out of range (#79 PR A)

    Code
      cat(msgs$fetwfe)
    Output
      q > 0 is not TRUE

