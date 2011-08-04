import bs

hdm = bs.bs_bos_core_data_storage.hdm()
hdm.get_keyword_manager().init(hdm)

hdm.read_keyword_file('D://projects//tests//test-models//Example-01//MODEL.DATA')
