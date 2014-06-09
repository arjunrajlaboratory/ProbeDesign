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


class FindProbesInput(ClassSerializer):
    __namespace__ = "findprobesinput"
    probeprefix = String(min_occurs=1,nillable=False)
    inseq = String(min_occurs=1,nillable=False)
    noligos = Integer(min_occurs=1,nillable=False)
    spacer_length = Integer(min_occurs=1,nillable=False)
    oligo_length = Integer(min_occurs=1,nillable=False)
    maskingflag = Boolean(min_occurs=1,nillable=False)
    species = String(min_occurs=1,nillable=False)

class ProbeOligo(ClassSerializer):
    __namespace__ = "probeoligo"
    GC = Decimal(min_occurs=1,nillable=False)    # percent GC content
    outseq = String(min_occurs=1,nillable=False) # 'ACTGTGTCATAGCT'
    label = String(min_occurs=1,nillable=False)  # 'HOTAIR_1'

class Alignment(ClassSerializer):
    __namespace__ = "alignment"
    raw_sequence = String(min_occurs=1,nillable=False)  # 'ACTGTGTCATAGCT'
    probe_oligos = String(min_occurs=1,nillable=False)  # '   CACAGTAT   '
    labels       = String(min_occurs=1,nillable=False)  # '   HOTAIR_1   '
    
class FindProbesOutput(ClassSerializer):
    __namespace__ = "findprobesoutput"
    oligos = Array(ProbeOligo)
    alignment = Alignment

class FISHProbeDesigner(DefinitionBase):
    """Set of SOAP web service methods to access probe design software."""
    
    @rpc(FindProbesInput,_returns=FindProbesOutput)
    def design_probes(self,input):
        # setup some log formating so we have a record of activity
        time_str = strftime("%a, %d %b %Y %H:%M:%S", localtime())
        print "Design probes request: %s" % time_str

        inseq = fasta.Fasta(input.inseq,strflag=True)
        
        results = find_probes.design(inseq,
                                    input.noligos,  
                                    input.oligo_length,
                                    input.spacer_length,
                                    input.maskingflag,
                                    input.species)

        alignment = probe_design.alignOutput(results['masked_seq'],
                                     results['output'][-1][1],input.oligo_length)
        # We have to replace ending spaces with non-whitespace characters
        # so the XML serializing wont cut out all the work we did to 
        # nicely align the positions of the probes and labels against inseq
        align = Alignment() 
        align.raw_sequence = '/' + inseq.one_line() + '/'  # masked is alignment[0]
        align.probe_oligos = '/' + alignment[1][0:-1] + '/'
        align.labels = '/' + alignment[2][0:-1] + '/'

        maxoligos = results['output'][-1] # scores, [1] matches [2] oligos
        maxoligos = probe_design.probeNames(maxoligos[2],input.probeprefix)
        olis = []
        for i in range(0,len(maxoligos)):  # for each designed oligo
            oli = ProbeOligo()
            oli.GC = maxoligos[i][0]
            oli.outseq = maxoligos[i][1]
            oli.label = maxoligos[i][2]
            olis.append(oli)

        output = FindProbesOutput()
        output.alignment = align
        output.oligos = olis
        
        print "Successful design of %d oligos" % len(olis)
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
        design_app = ValidatingApplication([FISHProbeDesigner],'tns')
        from cherrypy.wsgiserver import CherryPyWSGIServer
        server = CherryPyWSGIServer(('158.130.146.251',8080),
                                    design_app,
                                    server_name='localhost')
        server.start()
    except KeyboardInterrupt:
        server.stop()
    
