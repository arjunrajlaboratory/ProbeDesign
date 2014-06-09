"""SOAP server to host FISH probe design services.

Description:
    In order to provide access to FISH probe design software over the web,
    this program hosts a SOAP-WSDL web service accessed using client software
    that communicates via SOAP XML messages.

    To get the description of available methods and in/output formatting, use:
        http://hostname:port/?wsdl
    and the server will return the WSDL XML document. WSDL (web service
    description language) documents are like contracts that spell out the details
    of how request and response messages are formatted. With a valid WSDL, 
    client software can be have functions and data-types created "automagically" 
    by a wsdl2* tool where * is your favorite programming language.

Authors:
    Marshall J. Levesque 2011
    Arjun Raj 2011

"""

# modules specific to FISH probe design
import find_probes
import fasta
import bowtie_search
import probe_design

# modules for the SOAP server
from soaplib.service import rpc
from soaplib.service import DefinitionBase
from soaplib.serializers.primitive import String, Integer, Boolean, Decimal

from soaplib.wsgi import ValidatingApplication
from soaplib.serializers.clazz import ClassSerializer,Array
from soaplib.serializers import exception

# general modules for extra stuff
from time import localtime, strftime
import sys


class BowtieScreenInput(ClassSerializer):
    __namespace__ = "bowtiescreeninput"
    inseq = String(min_occurs=1,nillable=False)
    mer_length = Integer(min_occurs=1,nillable=False)
    seq_database = String(min_occurs=1,nillable=False)
    
class BowtieScreenOutput(ClassSerializer):
    __namespace__ = "bowtiescreenoutput"
    hits = Array(Integer)

class BowtieScreener(DefinitionBase):
    """Set of SOAP web service methods to access probe design software."""
    
    @rpc(BowtieScreenInput,_returns=BowtieScreenOutput)
    def screen_seqence(self,input):
        # setup some log formatting so we have a record of activity
        time_str = strftime("%a, %d %b %Y %H:%M:%S", localtime())
        print "Screen sequence request: %s" % time_str

        inseq = fasta.Fasta(input.inseq,strflag=True)
        '''
        results = find_probes.design(inseq,
                                    input.noligos,  
                                    input.oligo_length,
                                    input.spacer_length,
                                    input.maskingflag,
                                    input.species)
        '''

        '''
        output = FindProbesOutput()
        output.alignment = align
        output.oligos = olis
        '''
        
        output = BowtieScreenOutput()
        hts = bowtie_search.align_for_hits(inseq,input.mer_length,input.seq_database)
        output.hits = hts
        
        print "Successful screen"
        return output

    """ !!! THIS SECTION OF CODE USES SOAPLIB HOOKS API, but IT DOESN'T WORK !!! 

    def on_call(self,environ):
        # setup some log formating so we have a record of activity
        time_str = strftime("%a, %d %b %Y %H:%M:%S", localtime())
        msg = "Design probes request: %s" % time_str
        request.additional['tstart'] = msg
    
    def on_method_exception_object(self,environ,exc):
        msg = "Problem with method execution: %s" % repr(exc)
        request.additional['error'] = msg

    def on_return(self,environ,returnString):
        tstart = request.additional['tstart']
        if 'error' in request:
            print request.additional['error']
        else:
            print "Successful"
    """
    
if __name__=='__main__':
    try:
        """
        from wsgiref.simple_server import make_server
        server = make_server('158.130.146.251', 8080, Application([FISHProbeDesigner], 'tns'))
        server.serve_forever()
        """
        design_app = ValidatingApplication([BowtieScreener],'tns')
        from cherrypy.wsgiserver import CherryPyWSGIServer
        #server = CherryPyWSGIServer(('158.130.146.251',8080),
        server = CherryPyWSGIServer(('158.130.146.251',8000),
                                    design_app,
                                    server_name='localhost')
        server.start()
    except KeyboardInterrupt:
        server.stop()
    
