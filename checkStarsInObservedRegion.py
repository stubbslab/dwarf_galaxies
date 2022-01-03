from ObservedGalaxyStarData import ObservedGalaxyStarData
from DwarfGalDataArchive import DwarfGalDataArchive

def checkStarsInObservedRegion(population, pop_selection_method = 'none',fields_of_interest = [],masking_fields = []):
    dsph_archive = DwarfGalDataArchive()
    star_data = ObservedGalaxyStarData(population, pop_selection_method = 'none') 
    fields = star_data.Field
    corr_ra_one_pop = star_data.corrRa
    corr_dec_one_pop = star_data.corrDec
    if len(fields_of_interest) > 0:
        corr_ra_one_pop_untrimmed = corr_ra_one_pop[:]
        corr_dec_one_pop_untrimmed = corr_dec_one_pop[:]
        fields_untrimmed = fields
        corr_ra_one_pop = []
        corr_dec_one_pop = []
        fields = []
        for i in range(len(fields_untrimmed)):
            field = fields_untrimmed[i]
            if field in fields_of_interest:
                corr_ra_one_pop = corr_ra_one_pop + [corr_ra_one_pop_untrimmed[i]]
                corr_dec_one_pop = corr_dec_one_pop + [corr_dec_one_pop_untrimmed[i]]
                fields = fields + [fields_untrimmed[i]]

    obsLog = dsph_archive.obsLog[population[0]] 
    for i in range(len(corr_ra_one_pop)):
        inRegion = 0
        for log in obsLog:
            if len(masking_fields) == 0 or log['name'] in masking_fields:
                inRegion = inRegion + dsph_archive.inObservation(corr_ra_one_pop[i],corr_dec_one_pop[i],log)
        if inRegion  == 0:
            print 'Element ' + str(i) + ' located at (' + str(corr_ra_one_pop[i]) + ', ' + str(corr_dec_one_pop[i]) + ') not in region.  It is in field ' + fields[i] + '.'
