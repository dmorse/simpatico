RESTART             16000
FINISH

SET_CONFIG_IO  McConfigIo
READ_CONFIG        config
THERMALIZE            1.0
SET_CONFIG_IO  MdConfigIo
WRITE_CONFIG       config
WRITE_PARAM         param

