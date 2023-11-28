rule check_sv:
    output:
        sv=config["file"]["valid_sv"],
        repeat=config["file"]["valid_repeat"],
    run:
        get_valid_sv(insert_df, output.sv)
        get_valid_repeatmasker(insert_df, output.repeat)

